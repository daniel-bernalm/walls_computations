import numpy as np
from fractions import Fraction as fr
import math
from math import ceil, floor
import time, datetime
from multiprocessing import Pool

def in_ZZ(Fraction):
    return Fraction.denominator==1

def verify(tuple,beta,amb_sp):
    if amb_sp=="P3":
        r = tuple[0]
        c = tuple[1]
        d = tuple[2]
        e = tuple[3]
        a = d - fr(c**2, 2) + (beta)*c*(1-r) + (fr(beta**2, 2))*r*(1-r)
        b = e - fr(c, 6) + (beta)*(d - fr(r, 6)) + (fr(beta**2, 2))*c + (fr(beta**3, 6))*r
        c = 2*e - c*d + fr(c**3, 6) + (beta)*(d*(2-r) + (c**2)*(3*r-1)) + fr(beta**2, 2)*c*(2 + r*(r-3)) + fr(beta**3, 6)*r*(r-1)*(r-2)
    elif amb_sp=="Ab3":
        r = tuple[0]
        c = tuple[1]
        d = tuple[2]
        e = tuple[3]
        a = r*(beta**3) + 3*(beta**2)*c + 3*beta*d + e
        b = r*(beta**2) + 2*beta*c + d
        c = r*beta + c
    return in_ZZ(a) and in_ZZ(b) and in_ZZ(c) 

def grandverify(list,beta, amb_sp):
    a = len(list)
    granddefinitive = []
    #walls = []
    if amb_sp=="P3":
        for i in range(a):
            r1 = math.ceil(list[i][0])
            r2 = math.floor(list[i][1])
            definitive = ([], list[i][2],list[i][3],list[i][4])
            wl =  math.sqrt(fr(6*list[i][4], list[i][2])) 
            # This if comes from a conjecture. Change to 1 or 0 if needed. 
            if wl > (1):
                for j in range(r1,r2+1):
                    c=(j,list[i][2],list[i][3],list[i][4])
                    if verify(c,beta, amb_sp):
                        definitive[0].append(j)
            if len(definitive[0])>0:
                #walls.append(wl)
                granddefinitive.append(definitive)
    # probably not needed (CHECK)
    elif amb_sp=="Ab3":
        for i in range(a):
            r1 = math.ceil(list[i][0])
            r2 = math.floor(list[i][1])
            definitive = ([], list[i][2],list[i][3],list[i][4])
            wl = fr(list[i][4], list[i][2]) 
            for j in range(r1, r2+1):
                    c=(j,list[i][2],list[i][3],list[i][4])
                    if verify(c, beta, amb_sp):
                        definitive[0].append(j)
            if len(definitive[0])>0:
                #walls.append(wl)
                granddefinitive.append(definitive)
    return granddefinitive

class Sheaf:
    def __init__(self, R, D, k, amb_sp):
        self.R = -R
        self.D = D
        self.chern = (R, 0, D, 0)
        self.k = k
        self.amb_sp = amb_sp
        if k==1:
            self.beta = 0
        else:
            self.beta = fr(1, k)
        self.walls = set()
        self.dest = []
        self.maxwall = 0
        self.maxobject = []
    
    def possible_c6e(self):
        if self.amb_sp=="P3":
            solutions = []
            for d in range(1, 2*((self.k)**2)*self.D):
                trued = fr(d, 2*((self.k)**2))
                limit1 = 4*(trued**2)
                limit2 = 4*((self.D - trued)**2)
                min_limit = math.floor( ((self.k)**4)*min(limit1,limit2) )

                for c in range(1, min_limit + 1):
                    for e in range(1, min_limit + 1):
                        if 0 < c * e <= min_limit:
                            truec = fr(c, self.k)
                            truee = fr(e, 6*((self.k)**3))
                            # Ranks
                            rank1 = fr(-truec*(2*self.D - 2*trued), 6*truee) - self.R
                            rank2 = fr(truec*2*trued,6*truee)
                            solutions.append((rank1,rank2,truec,trued,truee))
        elif self.amb_sp=="Ab3": 
            solutions = []
            for d in range(1, ((self.k)**2)*self.D):
                trued = fr(d, ((self.k)**2))
                limit1 = trued**2
                limit2 = (self.D - trued)**2
                min_limit = math.floor( ((self.k)**4)*min(limit1,limit2) )

                for c in range(1, min_limit + 1):
                    for e in range(1, min_limit + 1):
                        if 0 < c * e <= min_limit:
                            truec = fr(c, self.k)
                            truee = fr(e, (self.k)**3)
                            # Ranks
                            rank1 = fr(-truec*(self.D - trued), 6*truee) - self.R
                            rank2 = fr(truec*trued,6*truee)
                            solutions.append((rank1,rank2,truec,trued,truee))
        return solutions

    def num_dest(self, number_division):
        list = self.possible_c6e()
        a = len(list)
        lenlist = math.floor( fr(a, number_division) )
        listblocks = [list[x:x+lenlist] for x in range(0, a, lenlist)]
        b = len(listblocks)

        if __name__ == "__main__":
            with Pool(processes=8) as pool:
                results = [pool.apply_async(grandverify, (listblocks[i], self.beta, self.amb_sp)) for i in range(b)]
                output = [res.get() for res in results]
                # Non-optimized but usually don't take a lot of time
                output2 = [x for x in output if x!=[]]
                self.dest = [output2[i][j] for i in range(len(output2)) for j in range(len(output2[i]))]
        return self.dest

t1 = time.time()
a = Sheaf(0,4,1, "Ab3")
print(a.num_dest(1))
t2 = time.time()
total = t2-t1
print(str(datetime.timedelta(seconds=total)))