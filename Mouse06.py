# Mouse T-cell repertoire model.   GL and MC 2016
# Thymic production plus death and division in the periphery.
# Age-dependence of thymus and body mass from Hogan et al. PNAS (2015)
# Time is measured in days.
# Gillespie algorithm
# Use scale factor 0.2 or less, if you are not on a machine with lots of RAM
from __future__ import print_function  # use python3 notation print()
from random import randrange,random,choice
from numpy import log,searchsorted,array
from math import exp
import matplotlib
matplotlib.use('Agg')
import pylab,copy,datetime

''' Model of the T-cells in a mouse '''

class Cell(object):
   ''' setup so that cell types inherit a counter '''
   number = 0
   def __init__(self):
      type(self).number += 1
   def __del__(self):
      type(self).number -= 1

class T(Cell):
   ''' T cell class '''
   def __init__(self):
      Cell.__init__(self)

class CD4(T):
   ''' CD4 T cell class '''
   def __init__(self):
      T.__init__(self)

class CD8(T):
   ''' CD8 T cell class '''
   def __init__(self):
      T.__init__(self)

cellnames=['CD4','CD8']
celltypes={'CD4':CD4,'CD8':CD8}

def death(thistype):
   ''' a cell dies '''
   thislist = celllists[thistype]
   thiscell = thislist.pop(randrange(len(thislist)))
   del thiscell

def death2(thistype):
   ''' a cell dies '''
   thiscell = celllists[thistype][-1]
   del celllists[thistype][-1]
   del thiscell

def division(thistype):
   ''' a cell divides '''
   thislist = celllists[thistype]
   thislist.append(celltypes[thistype]())

def thymus(thistype):
   ''' of a new clonotype, with nthy cells '''
   thislist = celllists[thistype]
   thislist.extend([celltypes[thistype]() for i in range(nthy)])

events={0:death2, 1:division, 2:thymus}

def makerates():
   ''' construct a list of rates of the Gillespie step '''
   rates = []
   ntotal = CD4.number + CD8.number
   aux=0
   for celltype in celltypes:
      aux = aux + mu[celltype]*celltypes[celltype].number
      rates.append(aux)
      aux = aux + gamma*celltypes[celltype].number/ntotal
      rates.append(aux)
      aux = aux + theta[celltype]
      rates.append(aux)  
   return rates,aux

def printinfo():
   print('At',int(t/7),'weeks, there are',CD4.number,'CD4 cells and',CD8.number,'CD8 cells.')
   print(' Daily thymic production is',int(thyout),'and peripheral division',int(gamma))


################## parameter values ####################################
# scaling factor
sfac = 0.10  # fraction of the whole mouse
# thymic production:
nthy,f8,thymax = 4,1.0/3,1000000*sfac  # cells per new clone, fraction CD8, daily thymic rate at 8 weeks
nu = 0.004 # decay rate
# periphery:
gmax = 200000*sfac
mu = {'CD4':0.030,'CD8':0.015}  # death rate
########################################################################
ncells = int(10*thymax) # initial cell number
CD4list = [CD4() for i in range(ncells//2)]
CD8list = [CD8() for i in range(ncells//2)]
celllists={'CD4':CD4list,'CD8':CD8list}
t,tmax = 0.0,63*7.
tforplot,cd4forplot,cd8forplot = [t],[CD4.number],[CD8.number]
gamma,thyout = 0,thymax*exp(-nu*(t-56))
print('scale factor =',sfac)
printinfo()

while t < tmax:
   gamma = gmax*(1-exp(-10*nu*t)) # peripheral division
   thyout = thymax*exp(-nu*(t-56))
   theta = {'CD4':thyout*(1-f8)/nthy,'CD8':thyout*f8/nthy}  # thymic production
   rates,sumrates = makerates()
   i = searchsorted(rates,random()*sumrates)
   events.get(i%3)(cellnames[i//3])
   t -= log(random())/sumrates
   if int(t*24) != int(tforplot[-1]*24):   
      tforplot.append(t)
      cd4forplot.append(CD4.number/sfac)
      cd8forplot.append(CD8.number/sfac)
      if int(tforplot[-1]/7) != int(tforplot[-2]/7):
         printinfo()

pylab.step(tforplot,cd4forplot,label='CD4 cells')
pylab.step(tforplot,cd8forplot,label='CD8 cells')
#pylab.step(tforplot,array(cd8forplot)+array(cd4forplot),label='total')
pylab.plot(tforplot,[7*thymax*exp(-nu*(t-56))/sfac for t in tforplot],':k',label='weekly thymic')
pylab.plot(tforplot,[7*gmax*(1-exp(-10*nu*t))/sfac for t in tforplot],':r',label='weekly division')
pylab.xlabel('$t$ / weeks')
pylab.xlim([0,max(tforplot)])
pylab.ylim([0,max(cd4forplot)*1.1])
pylab.xticks([28*i for i in range(int(tmax/28)+1)],[4*i for i in range(int(tmax/14)+1)])
pylab.yticks([1e7,2e7],[10,20])
pylab.ylabel('millions of cells')
pylab.legend(loc='upper right')
mystring = '$\\times$ '+str(sfac)+'   $\\nu=$'+str(nu)
pylab.title('Mouse'+mystring)
pylab.rcParams.update({'font.family': 'serif'})
todayraw = datetime.datetime.today()
mydate = todayraw.strftime("%d%b%H%M")
pylab.savefig('Mouse'+str(int(100*sfac))+str(mydate)+'.png')
pylab.show()

####################################################
# use // for integer division, works in python2 and python3
# this version does not work under python3
