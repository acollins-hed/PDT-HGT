import math
import numpy as np
import random
import time
from sets import Set
import argparse

parser = argparse.ArgumentParser(description="It's HGT")
parser.add_argument("--hgt", type=bool, default=True, help="It's HGT (default: 1).")
parser.add_argument("-t", "--time_steps", type=int, default=10, help="It's the time steps, must be an integer (default: 10).")
parser.add_argument("-p", "--pop_size", type=int, default=10, help="It's population size, must be an integer (default: 10).")
parser.add_argument("-nc", "--codon_number", type=int, default=8, help="It's codon number, must be an integer (default: 8).")
parser.add_argument("-na", "--aa_number", type=int, default=4, help="It's amino acid number, must be an integer (default: 4).")
parser.add_argument("-m", "--mu", type=float, default=0.0001, help="It's the mutation rate mu, must be in [0,1] (default: 0.0001).")
parser.add_argument("-n", "--nu", type=float, default=0.01, help="It's the mistranslation rate nu, must be in [0,1] (default: 0.01).")
parser.add_argument("-f", "--phi", type=float, default=0.99, help="It's the chemical distance scale between amino acids, phi, must be in (0,1) (default: 0.99).")
parser.add_argument("-H", "--mosaic", type=float, default=0.4, help="It's the fraction of the acceptor genome that is a mosaic due to HGT, must be in [0,1] (default: 0.4).")
parser.add_argument("-k", "--compatibility", type=float, default=1., help="It's the compatibility between donor and acceptor codes, must be in [0,1] (default: 1.0).")
args = parser.parse_args()

class HGT(object):
    #def __init__(self, initial_code, site_types, mu, message_mutation_matrix, pop_size, rng):
    """ 
    steps is an integer, it's the number of time steps
    trnas is the number of tRNAs and thus the number of codons as well
    aars is the number of aaRSs and thus the number of aa and site types
    mu is the mutation rate
    nu is the mistranslation rate
    phi is a chemical parameter used to calculate fitness
    H is a measure of the how much of the genome is mosaic
    pop_size is the number of individuals in the population
    rng is not used yet
    hgt is a Boolean variable. If true HGT happens, if false it doesn't
    """
    def __init__(self, steps, trnas, aars, mu, nu, phi, H, pk, pop_size, rng, hgt):
        if pop_size < 1:
            raise ValueError("The population size must be greater than 0.")
        
        if not 0 <= mu <= 1:
            raise RuntimeError("Mu must be in the range [0,1].")
        if not 0 <= nu <= 1:
            raise RuntimeError("Nu must be in the range [0,1].")
        if not 0 <= pk <= 1:
            raise RuntimeError("Nu must be in the range [0,1].")

        if not 0 < phi < 1:
            raise RuntimeError("Phi must be in the range (0,1).")
        #self._fitness_matrix = site_types.fitness_matrix
        #self._message_mutation_matrix = message_mutation_matrix
        #self._site_weights = site_types.weights
        self._pop_size = pop_size
        self._rng = rng
        self._mu = mu
        self._trnas = trnas
        self._H = H
        self._pk = pk
        self._aars = aars
        self.codes = None
        self._nu = nu
        self._phi = phi
        self._bits = int(np.log2(self._trnas))
        self._steps = steps
        self._hgt = hgt
    """ code generator generates a genetic code. It makes sure each aa has
    a codon. I likely will be changing this soon since the model I've been
    using as a template does not do it this way."""
    def code_generator(self):
        self.codes = np.zeros([self._pop_size,self._trnas])
        """for i in xrange(self._pop_size):
            for j in xrange(self._trnas):
                if j < self._aars:
                    self.codes[i,j] = j
                else:
                    self.codes[i,j] = np.random.random_integers(0,self._aars-1)
            self.codes[i,] = np.random.permutation(self.codes[i,])"""
        for j in xrange(self._trnas):
            if j < self._aars:
                self.codes[0,j] = j
            else:
                self.codes[0,j] = np.random.random_integers(0,self._aars-1)
        self.codes[0,] = np.random.permutation(self.codes[0,])
        for i in xrange(self._pop_size):
            self.codes[i,] = self.codes[0,].copy()
        return self.codes
    """ It retrieves an aa for a certain codon."""
    def get_aa(self,n,c):
        if c in self.codes[n,]:
            i = np.where(self.codes[n,] == c)[0][0]
            return i
        else:
            return -1
    """ Produces a set of nearest neighbors for a given codon."""
    def nearest_neighbors(self,from_):
        a = set([from_ ^ (1 << i) for i in xrange(self._bits)])
        b = set([from_ ^ (3 << (2*i)) for i in xrange(int(math.floor(float(self._bits)/2.)))])
        return a.union(b)
    """ This function returns the Hamming distance between different
    codes but as of now it does not work the way it's intended.
    I will attempt to fix it soon."""
    def code_hamming_distance(self,n,codes,teh_trna,old_aa,new_aa):
        dif = 0
        for i in xrange(self._pop_size):
            if i != n:
                if codes[i,teh_trna] == old_aa:
                    dif += 1
                if codes[i,teh_trna] == new_aa:
                    dif -= 1
        return dif
    """ swap used to swap aa between neighboring codons but now it 
    randomly mutates so a codon codes for a different aa. First, it checks
    that this won't irradicate an aa from a genetic code. I will change
    this soon as I don't think this is the best way to do this."""
    def swap(self,code):
        a = np.random.random_integers(0,self._trnas-1)
        b = code[a].copy()
        code[a] = np.random.random_integers(0,self._aars-1)
        '''if len(np.where(code == b)[0]) > 1:
            #c = np.random.random_integers(0,self._bits-1)
            code[a] = code[self.nearest_neighbors(a).pop()]'''
        return [a,b,code[a]]
        #a = np.random.random_integers(0,self._trnas-1)
        #b = np.random.random_integers(0,self._bits-1)
        #f = self.nearest_neighbors(a)[b]
        #d = code[a]
        #code[a] = code[f]
        #code[f] = d
    """ Calculates fitness of genetic code by looking at the chemical
    similarity among aas among neighboring codons."""
    def fitness(self,n,L,T,W,u,code):
        f = 1
        _sum = 0
        s = 1
        for i in xrange(self._trnas):
            for j in xrange(self._aars):
                for m in xrange(self._trnas):
                    #_sum += T[i,m]*W[self.get_aa(n,m),j]
                    _sum += T[i,m]*W[int(code[m]),j]
                s *= _sum**(L[j]*u[n,j,i])
                _sum = 0
            f *= s
            s = 1
        return f
    def N_c(self,c,c1):
        if c1 in self.nearest_neighbors(c):
            return 1
        else:
            return 0
    """ Calculates aa chemical distance score among neighboring codons."""
    def aa_dist(self,S,codes):
        mus = np.zeros([self._pop_size])
        for k in xrange(self._pop_size):
            for i in xrange(self._trnas):
                for j in range(i,self._trnas):
                    mus[k] += self.N_c(i,j)*S[int(codes[k,i]),int(codes[k,j])]
        return mus
    
    """ This was the model I wrote originally but I scrapped it. Continue to the
    model "vets" for the actual model."""
    def vetsigian(self,codes,u,HGT):
        
        if self.codes is None:
            self.code_generator()
        print self.codes[0]
        da_codes = self.codes.copy()
        
        data = np.zeros([2*self._steps,self._pop_size+2])
        #print data
        codons = set(range(0,self._trnas))

        # u is pop size by aaRS by tRNAs
        """ u is the different codon usage for each aa"""
        '''u = np.zeros([self._pop_size,self._aars,self._trnas])
        for i in xrange(self._pop_size):
            for j in xrange(self._aars):
                u[i,j,] = np.random.dirichlet([1]*self._trnas)
        uh = u.copy()'''

        """ A is a vector of aa chemical distances."""
        A = np.arange(0,1,1/float(self._aars))
        #print A

        """ W is a matrix of fitnesses for aa for each site type."""
        W = np.zeros([self._aars,self._aars])
        for i in xrange(self._aars):
            for j in xrange(self._aars):
                W[i,j] = self._phi**abs(A[i]-A[j])
        #print W
        #ssum = 0
        """ S is the aa similarity matrix."""
        S = np.zeros([self._aars,self._aars])
        mus = 0
        for i in xrange(self._aars):
            for j in xrange(self._aars):
                for k in xrange(self._aars):
                    mus += abs(W[i,k] - W[j,k])
                S[i,j] = mus.copy()
                #ssum += S[i,j]
                mus = 0
        #print S
        #print ssum
        """ T is a matrix of probablilities that a codon will be read as a
        neighboring codon."""
        T = np.zeros([self._trnas,self._trnas])
        for i in xrange(self._trnas):
            for j in xrange(self._trnas):
                if i == j:
                    T[i,j] = 1-self._nu
                elif j in self.nearest_neighbors(i):
                    T[i,j] = self._nu/(float(self._bits)+math.floor(float(self._bits)/2.))
                else:
                    T[i,j] = 0
        #print T
        """ pl is a vector of site types and L is their relative frequencies."""
        pl = np.random.dirichlet([1]*self._aars)
        L = pl/max(pl)
        #print L
        """M = np.diag([1-2*self._mu]*self._trnas)
        for i in xrange(self._trnas):
            for j in xrange(self._trnas):
                if i == (j-1)%self._trnas or j==(i-1)%self._trnas:
                    M[i,j] = self._mu"""
        """ M is a mutation matrix."""
        '''M = np.diag([1-self._bits*self._mu]*self._trnas)
        for i in xrange(self._trnas):
            for j in xrange(self._trnas):
                if i in self.nearest_neighbors(j) or j in self.nearest_neighbors(i):
                    M[i,j] = self._mu'''
        M = np.zeros([self._trnas,self._trnas])
        for i in xrange(self._trnas):
            for j in xrange(self._trnas):
                if i == j:
                    M[i,j] = 1-self._mu
                elif j in self.nearest_neighbors(i):
                    M[i,j] = self._mu/(float(self._bits)+math.floor(float(self._bits)/2.))
                else:
                    M[i,j] = 0
        print M

        '''Beginning of no HGT run'''
        print "Beginning of no HGT simulation."
        Data = open("test_data.txt","w")
        for it in xrange(self._steps):
            n = 0
            k = 0
            while(n == k):
                n = np.random.random_integers(0,self._pop_size-1)
                k = np.random.random_integers(0,self._pop_size-1)

            '''if(self._hgt):
                u[n,] = ((1-self._H)*u[n,] + self._H*u[k,]).copy()'''
                        
            f0 = self.fitness(n,L,T,W,u,self.codes[n,])
            #print "f0 ",f0
            code2 = self.codes[n,].copy()
            code3 = code2.copy()

            '''Let me try something'''
            """This block of code searches the set of codons. First it 
            randomly selects a codon from the set of codons c. It checks
            that it's not the only codon that codes for its amino acid. If
            it doesn't, then the codon is removed from c. If it does, then
            the amino acid it codes for is randomly changed and then the
            codon is removed from c."""
            c = codons.copy()
            while(len(c)>0):
                b = random.sample(c,1)[0]
                #if len(np.where(code2 == code2[b])[0]) > 1:
                code2[b] = np.random.random_integers(0,self._aars-1)
                c.remove(b)
                f1 = self.fitness(n,L,T,W,u,code2)
                #print "f1 ",f1
                if f1 <= f0:
                    code2 = code3.copy()
      
                else:
                    break
            '''I tried something'''
            
            '''for i in xrange(self._trnas*10):
                teh_trna, old_aa, new_aa = self.swap(code2)
                f1 = self.fitness(n,L,T,W,u,code2)
                #print "f1 ",f1
                if f1 < f0:
                    code2 = code3.copy()
      
                else:
                    break'''
            self.codes[n,] = code2.copy()

            for j in xrange(self._aars):
                #for i in xrange(self._trnas):
                #    print self.codes[n,i] 
                F = np.diag([W[int(self.codes[n,i]),j] for i in xrange(self._trnas)])
                Q = np.dot(M,F)
                evals, evecs = np.linalg.eig(Q)
                u[n,j,] = abs(evecs[:,np.where(evals==max(evals))[0][0]])

            #if it == 0:
            ham = 0.
            for i in xrange(self._trnas):
                for j in xrange(self._pop_size):
                    for l in range(j,self._pop_size):
                        if l != j and self.codes[j,i] != self.codes[l,i]:
                            ham += 1.
            HAM = ham/((float(self._pop_size)*float(self._pop_size-1))/2.)
            #else:
            #    ham += self.code_hamming_distance(n,teh_trna,old_aa,new_aa)
            #    HAM = ham/((float(self._pop_size)*float(self._pop_size-1))/2.)
            score = self.aa_dist(S)
            #print score
            data[it,0] = it
            data[it,1:self._pop_size+1] = score
            data[it,self._pop_size+1] = HAM
            Data.write(str(data[it,:])+"\n")
            if(it%1000 == 0):
                print it," steps"
        Data.close()
        print self._steps," steps"

        '''Beginning of HGT run'''
        print "Beginning of HGT simulation."
        self.codes = da_codes.copy()
        print self.codes[0,]
        Data1 = open("test_data1.txt","w")
        for it in xrange(self._steps):
            n = 0
            k = 0
            while(n == k):
                n = np.random.random_integers(0,self._pop_size-1)
                k = np.random.random_integers(0,self._pop_size-1)

            uh[n,] = ((1-self._H)*uh[n,] + self._H*uh[k,]).copy()
            
            f0 = self.fitness(n,L,T,W,uh,self.codes[n,])
            #print "f0 ",f0
            code2 = self.codes[n,].copy()
            code3 = code2.copy()

            '''Let me try something'''
            """This block of code searches the set of codons. First it 
            randomly selects a codon from the set of codons c. It checks
            that it's not the only codon that codes for its amino acid. If
            it doesn't, then the codon is removed from c. If it does, then
            the amino acid it codes for is randomly changed and then the
            codon is removed from c."""
            c = codons.copy()
            while(len(c)>0):
                b = random.sample(c,1)[0]
                #if len(np.where(code2 == code2[b])[0]) > 1:
                code2[b] = np.random.random_integers(0,self._aars-1)
                c.remove(b)
                f1 = self.fitness(n,L,T,W,uh,code2)
                #print "f1 ",f1
                if f1 <= f0:
                    code2 = code3.copy()
      
                else:
                    break
            '''I tried something'''
            
            '''for i in xrange(self._trnas*10):
                teh_trna, old_aa, new_aa = self.swap(code2)
                f1 = self.fitness(n,L,T,W,u,code2)
                #print "f1 ",f1
                if f1 < f0:
                    code2 = code3.copy()
      
                else:
                    break'''
            self.codes[n,] = code2.copy()

            for j in xrange(self._aars):
                #for i in xrange(self._trnas):
                #    print self.codes[n,i] 
                F = np.diag([W[int(self.codes[n,i]),j] for i in xrange(self._trnas)])
                Q = np.dot(M,F)
                evals, evecs = np.linalg.eig(Q)
                uh[n,j,] = abs(evecs[:,np.where(evals==max(evals))[0][0]])

            #if it == 0:
            ham = 0.
            for i in xrange(self._trnas):
                for j in xrange(self._pop_size):
                    for l in range(j,self._pop_size):
                        if l != j and self.codes[j,i] != self.codes[l,i]:
                            ham += 1.
            HAM = ham/((float(self._pop_size)*float(self._pop_size-1))/2.)
            #else:
            #    ham += self.code_hamming_distance(n,teh_trna,old_aa,new_aa)
            #    HAM = ham/((float(self._pop_size)*float(self._pop_size-1))/2.)
            score = self.aa_dist(S)
            #print score
            data[it,0] = it
            data[it,1:self._pop_size+1] = score
            data[it,self._pop_size+1] = HAM
            Data1.write(str(data[it,:])+"\n")
            if(it%1000 == 0):
                print it," steps"
            
        Data1.close()
        print self._steps," steps"
        #print self.codes

    def vets(self,code,u,pl,HGT,fname):
        #print fname

        codes = np.zeros([self._pop_size,self._trnas])

        for i in xrange(self._pop_size):
            codes[i,] = code.copy()
        
        data = np.zeros([self._steps,self._pop_size+2])

        codons = set(range(0,self._trnas))
        
        
        #print u

        # u is pop size by aaRS by tRNAs

        """ A is a vector of aa chemical distances."""
        A = np.arange(0,1,1/float(self._aars))
        #print A

        """ W is a matrix of fitnesses for aa for each site type."""
        W = np.zeros([self._aars,self._aars])
        for i in xrange(self._aars):
            for j in xrange(self._aars):
                W[i,j] = self._phi**abs(A[i]-A[j])
        #print W

        """ S is the aa similarity matrix."""
        S = np.zeros([self._aars,self._aars])
        mus = 0
        for i in xrange(self._aars):
            for j in xrange(self._aars):
                for k in xrange(self._aars):
                    mus += abs(W[i,k] - W[j,k])
                S[i,j] = mus.copy()
                mus = 0
        #print S

        """ T is a matrix of probablilities that a codon will be read as a
        neighboring codon."""
        T = np.zeros([self._trnas,self._trnas])
        for i in xrange(self._trnas):
            for j in xrange(self._trnas):
                if i == j:
                    T[i,j] = 1-self._nu
                elif j in self.nearest_neighbors(i):
                    T[i,j] = self._nu/(float(self._bits)+math.floor(float(self._bits)/2.))
                else:
                    T[i,j] = 0
        #print T

        """ pl is a vector of site types and L is their relative frequencies."""
        L = (pl/max(pl)).copy()
        #print L

        M = np.zeros([self._trnas,self._trnas])
        for i in xrange(self._trnas):
            for j in xrange(self._trnas):
                if i == j:
                    M[i,j] = 1-self._mu
                elif j in self.nearest_neighbors(i):
                    M[i,j] = self._mu/(float(self._bits)+math.floor(float(self._bits)/2.))
                else:
                    M[i,j] = 0
        #print M

        #print fname
        Data = open(fname,"w")
        for it in xrange(self._steps):
            n = 0
            k = 0
            while(n == k):
                n = np.random.random_integers(0,self._pop_size-1)
                k = np.random.random_integers(0,self._pop_size-1)

            if(HGT):
                u[n,] = ((1-self._H)*u[n,] + self._H*u[k,]).copy()
                        
            f0 = self.fitness(n,L,T,W,u,codes[n,])
            #print "f0 ",f0
            code2 = codes[n,].copy()

            '''Let me try something'''
            """This block of code searches the set of codons. First it 
            randomly selects a codon from the set of codons c. It checks
            that it's not the only codon that codes for its amino acid. If
            it doesn't, then the codon is removed from c. If it does, then
            the amino acid it codes for is randomly changed and then the
            codon is removed from c."""
            c = codons.copy()
            while(len(c)>0):
                teh_rna = random.sample(c,1)[0]
                #if len(np.where(code2 == code2[b])[0]) > 1:
                code2[teh_rna] = np.random.random_integers(0,self._aars-1)
                c.remove(teh_rna)
                f1 = self.fitness(n,L,T,W,u,code2)
                #print "f1 ",f1
                if f1 <= f0:
                    code2 = codes[n,].copy()
      
                else:
                    break
            '''I tried something'''
            
            old_aa = codes[n,teh_rna].copy()
            new_aa = code2[teh_rna]
            codes[n,] = code2.copy()

            for j in xrange(self._aars):
                F = np.diag([W[int(codes[n,i]),j] for i in xrange(self._trnas)])
                Q = np.dot(M,F)
                evals, evecs = np.linalg.eig(Q)
                u[n,j,] = abs(evecs[:,np.where(evals==max(evals))[0][0]])

            if it == 0:
                ham = 0.
                for i in xrange(self._trnas):
                    for j in xrange(self._pop_size):
                        for l in range(j,self._pop_size):
                            if l != j and codes[j,i] != codes[l,i]:
                                ham += 1.
                #print ham
                HAM = ham/((float(self._pop_size)*float(self._pop_size-1))/2.)
            else:
                ham += self.code_hamming_distance(n,codes,teh_rna,old_aa,new_aa)
                '''ham1 = 0.
                for i in xrange(self._trnas):
                    for j in xrange(self._pop_size):
                        for l in range(j,self._pop_size):
                            if l != j and codes[j,i] != codes[l,i]:
                                ham1 += 1.
                #print ham - ham1
                HAM1 = ham1/((float(self._pop_size)*float(self._pop_size-1))/2.)'''
                HAM = ham/((float(self._pop_size)*float(self._pop_size-1))/2.)
                #print HAM - HAM1
            score = self.aa_dist(S,codes)
            #print score
            data[it,0] = it
            data[it,1:self._pop_size+1] = score
            data[it,self._pop_size+1] = HAM
            Data.write(str(data[it,:])+"\n")
            #print it, "steps"
            if(it%1000 == 0):
                print it," steps"
        Data.close()
        print self._steps," steps"
        #print "Final codes ",codes
        #print u
    def run(self):
        u = np.zeros([self._pop_size,self._aars,self._trnas])
        for i in xrange(self._pop_size):
            for j in xrange(self._aars):
                u[i,j,] = np.random.dirichlet([1]*self._trnas)
        U = u.copy()
        #print u
                
        code = np.zeros([self._trnas])
        for j in xrange(self._trnas):
            if j < self._aars:
                code[j] = j
            else:
                code[j] = np.random.random_integers(0,self._aars-1)
        code = np.random.permutation(code)
        #print "Starting code: ",code
        pl = np.random.dirichlet([1]*self._aars)
        print "Beginning of HGT simulation."
        self.vets(code,U,pl,True,"data_hgt.txt")
        #print u
        print "Beginning of no HGT simulation."
        self.vets(code,u,pl,False,"data_nohgt.txt")

sclock,scpu = [time.time(),time.clock()]
h = HGT(args.time_steps,args.codon_number,args.aa_number,args.mu,args.nu,args.phi,args.mosaic,args.compatibility,args.pop_size,3,args.hgt)
#h.code_generator()
#print h.codes,"\n"
#print h.codes[0,0],"\n"
h.run()
eclock,ecpu = [time.time(),time.clock()]
print "It took ",(eclock-sclock)/3600," clock hours and ",(ecpu-scpu)/3600," hours in cpu time.\n"
#print h.codes,"\n\n"
