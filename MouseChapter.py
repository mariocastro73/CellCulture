# Mouse T-cell repertoire model.   GL and MC 2016
# Thymic production plus death and division in the periphery.
# Age-dependence of thymus and body mass from Hogan et al. PNAS (2015)
# Time is measured in days. Gillespie algorithm
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
   def __init__(self,tcr=None):
      Cell.__init__(self)
      if tcr == None:
          self.tcr = randrange(0,1e12) # Possible different clonotypes
      else:
          self.tcr = tcr # If passed in the constructor, it's inherited

class CD4(T):
   ''' CD4 T cell class '''
   def __init__(self,tcr=None):
      T.__init__(self,tcr)

class CD8(T):
   ''' CD8 T cell class '''
   def __init__(self,tcr=None):
      T.__init__(self,tcr)


class CellPopulation(object):
   ''' Collection of cells and methods to manipulate them ''' 
   def __init__(self,ncells):
      self.cellnames=['CD4','CD8']
      self.celltypes={'CD4':CD4,'CD8':CD8}
      ###### The initial list of cells in the population #####
      self.CD4list = [CD4() for i in range(ncells//2)]
      self.CD8list = [CD8() for i in range(ncells//2)]
      self.celllists={'CD4':self.CD4list,'CD8':self.CD8list}
      ###### What methods do we have to manipulate cells #####
      self.events={0:self.death, 1:self.division, 2:self.thymus}

   def death(self,thistype):
      ''' a cell dies '''
      index = randrange(0,self.celltypes[thistype].number)
      thiscell = self.celllists[thistype][index]
      del self.celllists[thistype][index]
   
   def division(self,thistype):
      ''' a cell divides '''
      index = randrange(0,self.celltypes[thistype].number)
      tcr = self.celllists[thistype][index].tcr
      newcell = self.celltypes[thistype](tcr)
      self.celllists[thistype].append(newcell)
   
   def thymus(self,thistype):
      ''' creation of cells of a new clonotype, with nthy cells '''
      thislist = self.celllists[thistype]
      extendlist = [self.celltypes[thistype]() for i in range(nthy)]
      thislist.extend(extendlist)

class simulation(object):
   ''' Simulation class '''
   def __init__(self,tmax):
       self.tmax = tmax # Maximum simulation time
       self.t = 0 # Age of the mouse. Starting at 0
       self.tforplot = [self.t]# Auxiliary variables for output
       self.cd4forplot = [CD4.number]#Auxiliary variables for output
       self.cd8forplot = [CD8.number] #Auxiliary variables for output
       self.ratioforplot = [f8] #Auxiliary variables for output
       self.cellpop = CellPopulation(ncells) # Create ncells

   ''' Scheduler '''
   def scheduler(self):
       while self.t < self.tmax: # While not at tmax
          self.gamma = gmax*(1-exp(-10*nu*self.t)) # peripheral division
          self.thyout = thymax*exp(-nu*(self.t-56))# thymic production
          self.theta = {'CD4':self.thyout*(1-f8)/nthy,'CD8':self.thyout*f8/nthy}  
          rates,sumrates = self.makerates() # Create gillespie rates
          i = searchsorted(rates,random()*sumrates) # Find the chosen event
          self.cellpop.events.get(i%3)(self.cellpop.cellnames[i//3]) # Execute one event
          self.t -= log(random())/sumrates # Gillespie time update

          if int(self.t*24) != int(self.tforplot[-1]*24): # Store for output
             self.tforplot.append(self.t)
             self.cd4forplot.append(CD4.number/sfac)
             self.cd8forplot.append(CD8.number/sfac)
             self.ratioforplot.append(1.*CD4.number/CD8.number)
             if int(self.tforplot[-1]/7) != int(self.tforplot[-2]/7):
                print('At',int(self.t/7),'weeks, there are',CD4.number,
                        'CD4 cells and',CD8.number,'CD8 cells.')
                print(' Daily thymic production is',int(self.thyout),
                        'and peripheral division',int(self.gamma))

   def makerates(self):
      ''' construct a list of rates of the Gillespie step '''
      rates = []
      ntotal = CD4.number + CD8.number
      aux=0
      for celltype in self.cellpop.celltypes:
         aux = aux + mu[celltype]*self.cellpop.celltypes[celltype].number # Death
         rates.append(aux)
         aux = aux + self.gamma*self.cellpop.celltypes[celltype].number/ntotal # Division
         rates.append(aux)
         aux = aux + self.theta[celltype] # Create new clonotypes in the thymus
         rates.append(aux)  
      return rates,aux

   def visualization(self):
      ''' create a plot with the time course of CD4's and CD8's '''
      pylab.clf() # Clear the previous plot
      pylab.step(self.tforplot,self.cd4forplot,label='CD4 cells')
      pylab.step(self.tforplot,self.cd8forplot,label='CD8 cells')
      pylab.plot(self.tforplot,[7*thymax*exp(-nu*(t-56))/sfac for t in self.tforplot],
              ':k',label='weekly thymic')
      pylab.plot(self.tforplot,[7*gmax*(1-exp(-10*nu*t))/sfac for t in self.tforplot],
              ':r',label='weekly division')
      pylab.xlabel('$t$ / weeks')
      pylab.xticks([28*i for i in range(int(tmax/28)+1)],
                   [4*i for i in range(int(tmax/14)+1)])
      pylab.ylabel('millions of cells')
      pylab.yticks([1e7,2e7],[10,20])
      pylab.legend(loc='upper right')
      mydate = datetime.datetime.today().strftime("%d%b%H%M")
      pylab.savefig('Mouse-'+str(int(100*sfac))+str(mydate)+'.png')
      
      pylab.clf() # Clear the previous plot
      pylab.step(self.tforplot,self.ratioforplot,label='CD4/CD8 ratio')
      pylab.plot(self.tforplot,[1 for t in self.tforplot])
      pylab.xlabel('$t$ / weeks')
      pylab.xlim([0,max(self.tforplot)])
      pylab.xticks([28*i for i in range(int(tmax/28)+1)],
                   [4*i for i in range(int(tmax/14)+1)])
      pylab.ylim([0,max(self.ratioforplot)*1.1])
      pylab.ylabel('ratio CD4 to CD8')
      pylab.legend(loc='upper right')
      mydate = datetime.datetime.today().strftime("%d%b%H%M")
      pylab.savefig('Mouse-ratio'+str(int(100*sfac))+str(mydate)+'.png')

####################  global parameter values ####################################
# scaling factor
sfac = 0.1  # fraction of the whole mouse (use sfac = 1 for whole mouse)
# thymic production:
nthy,f8 = 4,1.0/5  # cells per new clone, fraction CD8
thymax,nu = 1000000*sfac,0.004 # daily thymic rate at 8 weeks, decay rate
# periphery:
gmax,mu = 200000*sfac,{'CD4':0.030,'CD8':0.015}  # peripheral division and death 
ncells = int(10*thymax) # initial cell number
tmax = 63*7 # Total number of days of the experiment
##################################################################################

sim = simulation(tmax) # Create a new simulation
sim.scheduler() # Start the scheduler
sim.visualization() # At the end of the simulation, produce the output

##################################################################################
# python2.7 is recommended
# use // for integer division, works in python2 and python3
##################################################################################
