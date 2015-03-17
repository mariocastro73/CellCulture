#############################################################
# March 2015
# T cell classes
#
##############################################################
from __future__ import print_function
from random import shuffle,sample,expovariate,random
from pylab import *
#import numpy as np
import copy

#ion()

class Cell(object):
    number = 0
    instances = []
    ''' store number and instances'''
    def __init__(self):
        type(self).instances.append(self)
        type(self).number += 1
        self.ct = ''
    def delete(self):
        type(self).instances.remove(self)
        type(self).number -= 1

class T(Cell):
   ''' T cell class. attributes are:
   cy cycling
   la TCR avidity
   ge generation
   mt mitosis time
   '''
   def __init__(self,thisg):
       Cell.__init__(self)
       self.cy = False
       self.la = 1.0
       if thisg==0: self.la=expovariate(24.0)
       self.ge = thisg
       self.mt = 13.0
       self.ct = 'T-cell'
   def __str__(self):
       return "%s // cy=%s la=%f ge=%f mt=%f"%(self.ct,self.cy,self.la,self.ge,self.mt)
   def entercycle(self):
       self.cy = True
       self.mt = 12.0
   def advance(self):
       ''' reset attributes just before division '''
       self.ge += 1
       self.cy = False
       self.mt = 13.0
   def celldivision(self):
       ''' use deepcopy, suggested by Mario'''
       self.advance()
       newcell = copy.deepcopy(self)
       type(self).instances.append(newcell)
       type(self).number += 1
       return newcell
   def cellstep(self,dt):
       if self.cy:
           self.mt -= dt
           if self.mt<0:
               return True
       elif stim*dt*self.la>random():
               self.entercycle()
       return False

class CD4(T):
    ''' CD4 T cell class '''
    def __init__(self,thisg):
        T.__init__(self,thisg)
        self.ct = 'CD4'

class CD8(T):
    ''' CD8 T cell class '''
    number = 0
    instances = []
    def __init__(self,thisg):
        T.__init__(self,thisg)
        self.ct = 'CD8'

celltype={'CD4':CD4,'CD8':CD8}

class CellPopulation(object):
    def __init__(self,ncells,ctype,initialcondition):
        self.ncells = ncells
        self.celllist = [celltype[ctype](initialcondition) for i in xrange(ncells)]    

    #@profile
    def step(self,dt):
        ''' in a time interval dt, some cells enter cycle and other complete it'''
        self.celllist += [cell.celldivision() for cell in self.celllist if cell.cellstep(dt)]

    def rsublist(self,n):
        '''random sublist of n cells, unweighted'''
        return set(sample(self.celllist,n))
        
    def isublist(self,n,wat):
        '''random sublist of n cells, weighted by integer attribute wat'''
        biglist = []
        for i in range(len(self.celllist)):
            biglist.extend([i]*getattr(self.celllist[i],wat))
        sublist = sorted(sample(biglist,n))
        sublist.reverse()
        return sublist
        
    def gennumbers(self,maxn):
        i=0
        g = []
        while i<maxn:
            g.append(len([thiscell.ge for thiscell in self.celllist if thiscell.ge==i]))
            i += 1
        return g


def myprint(mylist):
    [print("{:.2f}".format(x),end=' ') for x in mylist]
    print()

# now let's try it out
cdeights = CellPopulation(10000,'CD8',0)
stim = 1.0
t = 0
dt = 0.1
tmax =  72.0  # time is measured in hours
g = cdeights.gennumbers(8)
# set the graphics up
config,= step(arange(len(g)),g,where='post')
timetext = text(len(g)-2, max(g)*1.9,'t='+str(t))
ylim([0,2*max(g)])
xticks(arange(6)+0.5,arange(8))
xlabel('number of divisions')
ylabel('number of cells')
while t < tmax-dt/100:
    t += dt
    cdeights.step(dt)
    g = cdeights.gennumbers(8)
    #print(t,CD8.number,len(cdeights.celllist),len(CD8.instances),g)

    #config.set_ydata(g)
    #timetext.set_text(str(t)+' hours')
    #draw()

#bar(arange(len(g)),g)
#savefig('generations.png')

print(t,CD8.number,len(cdeights.celllist),len(CD8.instances),g)


