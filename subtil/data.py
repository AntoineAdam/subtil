import re
import traceback
import numpy as np
import matplotlib.pyplot as plt
from math import log

from .utils import exp_interpolation, remove_none, moving_average,is_float

_genealogy_regex = re.compile("^((descendantof)|(ancestorof)) +[0-9]+\.[0-9]+$")


#TODO: remove
_cell_arff_attribute_types={	'id':'numeric',
						'lineageid':'{}',
						'lefttreeposition':'numeric',
						'outtreeposition':'numeric',
						'index':'{}',
						'pole':'{0,1}',
						'generation':'numeric',
						'age':'numeric',
						'youth':'numeric',
						'birthtime':'numeric',
						'divisiontime':'numeric',
						'lifetime':'numeric',    #time in the experiment, not necessarily from brith to division
						'doublingtime':'numeric',
						'lagtime':'numeric',
						'length':'list',
						'instantlengthgrowthrate':'list',
						'avglength':'numeric',
						'lengthatbirth':'numeric',
						'lengthatdivision':'numeric',
						'lengthincrease':'numeric',
						'lengthgrowthrate':'numeric',
						'maxwidthgrowthrate':'numeric',
						'areagrowthrate':'numeric',
						'alive':'{0,1}',
						'timeatshock':'numeric',
						'lengthatshock':'numeric',
						'lengthincreaseatshock':'numeric',
						'lengthgrowthrateatshock':'numeric',
						'fluorescenceatshock':'numeric',
						'avgfluorescence':'numeric',
						'fluorescenceareaatshock':'numeric',
						'fluorescencefociatshock':'numeric',
						'cycleatshock':'numeric',
						'populationgrowthrate':'numeric',
						'areagrowthrate':'numeric',
						'heatshocktime':'numeric',
						'survivalrate':'numeric',
						'populationatshock':'numeric',
						'areaatshock':'numeric',
						'areagrowthrateatshock':'numeric',
						'survivalrate':'numeric'}
	
###################################################
# Classes                               																	             #
###################################################

# Cell class
# fluorescence is the signal from microbeTracker. multiplied by 2^16-1, it gives the total fluorescence for a cell at a given time. To get an average fluorescence of the cell, divide by the area at the same time.
class Cell(object):

	"""Cell object.

	Represents a cell object.

	Parameters
	-------------

	id: int
		unique id of the cell in the database
	lineageid: int
		id of the cell's lineage
	parent: Cell
		the parent Cell
	pole: int
		0 if the cell inherited the old pole of its parent, 1 otherwise.
	lefttreeposition: int
		position in the binary tree where oldpold cells are stored on the left.
	birthtime: double
		time in the experiment when the cell was born.
	divisiontime: double
		time in the experiment when the cell divided.
	time: array-like (double), optional
		timestamps in the experiment, in minutes.
	length: array-like (double), optional
		length of the cell at each timestamp, in pixels.
	maxwidth: array-like (double), optional
		maximum width of the cell at each timestamp, in pixels.
	area: array-like (double), optional
		area of the cell at each timestamp, in pixel-squares.
	fluorescence: array-like (double), optional
		fluorescence intensity of the whole cell at each timestamp, in ?.
	nucleoids: array-like (int), optional
		Number of fluorescence nucleoids.
	relativearea: array-like (double), optional
		Relative area of the cell covered by the nucleoids.
	
	"""

	attribute_types = {	'id':'nominal',
				'lineageid':'nominal',
				'lefttreeposition':'nominal',
				'outtreeposition':'nominal',
				'index':'treeindex',
				'pole':'nominal',
				'generation':'numeric',
				'age':'numeric',
				'youth':'numeric',
				'birthtime':'numeric',
				'divisiontime':'numeric',
				'lifetime':'numeric',
				'doublingtime':'numeric',
				'lagtime':'numeric',
				'length':'list',
				'instantlengthgrowthrate':'list',
				'avglength':'numeric',
				'lengthatbirth':'numeric',
				'lengthatdivision':'numeric',
				'lengthincrease':'numeric',
				'lengthgrowthrate':'numeric',
				'maxwidthgrowthrate':'numeric',
				'areagrowthrate':'numeric',
				'alive':'nominal',
				'timeatshock':'numeric',
				'lengthatshock':'numeric',
				'lengthincreaseatshock':'numeric',
				'lengthgrowthrateatshock':'numeric',
				'avgfluorescence':'numeric',
				'fluorescenceatshock':'numeric',
				'fluorescenceareaatshock':'numeric',
				'fluorescencefociatshock':'numeric',
				'cycleatshock':'numeric',
				'survivalrate':'numeric'}

	
	def __init__(self, id, lineageid, parent, pole, lefttreeposition, birthtime, divisiontime, time=np.array([]), length=np.array([]), maxwidth=np.array([]), area=np.array([]), fluorescence=np.array([]), nucleoids=np.array([]), relativearea=np.array([])):
		self.id=id
		self.lineageid=lineageid
		self.lineage=None
		self.parent=parent
		self.children=[]
		self.pole = pole
		self.lefttreeposition=lefttreeposition
		if lefttreeposition is None:
			self.generation=None
			self.age = None
			self.youth = None
		else:
			self.generation=int(log(lefttreeposition,2))
			age=0
			p=lefttreeposition
			while p%2==0 and p>3:
				age=age+1
				p=p/2
			self.age=age
			youth=0
			p=lefttreeposition
			while p%2==1 and p>3:
				youth=youth+1
				p=p/2
			self.youth=youth
		self.birthtime=birthtime
		self.divisiontime=divisiontime
		self.time=time		# in minutes
		self.length=length		# in pixels
		self.maxwidth=maxwidth
		self.area=area
		self.fluorescence=fluorescence
		self.nucleoids=nucleoids
		self.relativearea=relativearea

	
	def __str__(self):
		return "Cell {}: {}".format(self.id, self.get("index"))

	def get(self, attribute):
		
		s=attribute.split('_')
		
		try:
		
			if s[0]=='parent':
				return self._parentget(attribute[7:])
			elif s[0]=='sister':
				return self._sisterget(attribute[7:])
			elif s[0]=='daughter':
				return self._daughtersimilarity('_'.join(s[1:len(s)-1]), s[-1])
			elif hasattr(self, attribute):
				if attribute in self.__dict__:
					return getattr(self,attribute)
				else:
					return getattr(self,attribute)()
			else:
				print("attribute not found: {0}".format(attribute))
				return None
		
		except Exception as ex:
			print "Error while computing {} for cell {}".format(attribute, cell.id)
			print traceback.format_exc()
			return None


	def outtreeposition(self):
		if self.lefttreeposition<=3: #works if None
			return self.lefttreeposition
		else:
			parentposition = self.parent.get("outtreeposition")
			return parentposition*2 + ((parentposition % 2) if (self.pole==0) else (1-parentposition%2))

	def index(self):
		return "{0}.{1}".format(self.lineageid,self.lefttreeposition)

	def lifetime(self):
		return self.time[-1]-self.time[0]+self.lineage.timestep

	def doublingtime(self):
		if self.divisiontime is None or self.birthtime is None:
			return None
		else:
			return self.divisiontime-self.birthtime

	def lengthgrowthrate(self):
		return self._growthrate("length")
	
	def avglength(self):
		return np.mean(self.length)
	
	def lengthatbirth(self):
		if self.birthtime is None:
			return None
		return self.lengthat(self.time[0])

	def lengthatdivision(self):
		if self.divisiontime is None:
			return None
		return self.lengthat(self.time[-1])

	def lengthincrease(self):
		return self.lengthat(self.time[-1])-self.lengthat(self.time[0])

	def maxwidthgrowthrate(self):
		return self._growthrate("maxwidth")
	def areagrowthrate(self):
		return self._growthrate("area")

	def descendantalive(self):
		if hasattr(self,"alive"):
			return self.alive
		else:
			if len(self.children)==2:
				a1=self.children[0].get("descendantalive")
				a2=self.children[1].get("descendantalive")
				if a1 is not None and a2 is not None:
					return float(a1+a2)/2
				elif a1 is not None:
					return a1
				elif a2 is not None:
					return a2
				else:
					return None
			else:
				return None
	def timeatshock(self):
		if hasattr(self,"alive"):
			return self.time[-1]-self.time[0]
		else:
			return None
	def lengthatshock(self):
		if hasattr(self,"alive"):
			return self.length[-1]
		else:
			return None
	def lengthincreaseatshock(self):
		if hasattr(self,"alive"):
			return self.length[-1]-self.length[0]
		else:
			return None
	
	def lengthgrowthrateatshock(self):
		if hasattr(self,"alive"):
			return self.instantlengthgrowthrateat(self.time[-1])
		else:
			return None

	def avgfluorescence(self):
		if hasattr(self,"fluorescence"):
			ftime = np.array([t%self.lineage.fluotimestep==0 for t in self.time])
			return np.sum(np.array(self.fluorescence)[ftime]/np.array(self.area)[ftime])*65535/np.sum(ftime)
		else:
			return None
			
	def fluorescenceatshock(self):
		if hasattr(self,"alive"):
			return self.totalfluorescenceat(self.time[-1])
		else:
			return None
			
	def fluorescenceareaatshock(self):
		if hasattr(self,"alive"):
			return self.relativearea[-1]
		else:
			return None
	
	def fluorescencefociatshock(self):
		if hasattr(self,"alive"):
			return self.nucleoid[-1]
		else:
			return None
	
	def getsister(self):
		if self.parent is None:
			return None
		else:
			siblings=self.parent.children
			if len(siblings)==2:
				if self.id == siblings[0].id:
					return siblings[1]
				else:
					return siblings[0]
		return None
	
	def lengthat(self,t):
		if not t in self.time:
			return None
		return self.length[np.where(self.time==t)[0][0]]
	
	def areaat(self,t):
		if not t in self.time:
			return None
		return self.area[np.where(self.time==t)[0][0]]
	
	def totalfluorescenceat(self,t):
		if not t in self.time or t%self.lineage.fluotimestep != 0:
			return None
		return self.fluorescence[np.where(self.time==t)[0][0]]*65535 if self.fluorescence[np.where(self.time==t)[0][0]] is not None else None
	
	def avgfluorescenceat(self,t):
		if not t in self.time or t%self.lineage.fluotimestep != 0:
			return None
		i=np.where(self.time==t)[0][0]
		return self.fluorescence[i]*65535/self.area[i] if self.fluorescence[i] is not None else None
	
	# does an exponential interpolation of the length
	def _growthrate(self, attr='length'):
		if attr=='length':
			if hasattr(self, 'lengthgrowthrate_'):
				return self.lengthgrowthrate_
			values=self.length
		elif attr=='maxwidth':
			if hasattr(self, 'maxwidthgrowthrate_'):
				return self.maxwidthgrowthrate_
			values=self.maxwidth
		elif attr=='area':
			if hasattr(self, 'areagrowthrate_'):
				return self.areagrowthrate_
			values=self.area
		else:
			print("wrong attribute for growthrate: {0}".format(attr))
			return 0
		
		if len(values)<=1:
			print("warning: only one value for {0} of cell {1}".format(attr,self.id))
			self.lengthgrowthrate_=0
			self.maxwidthgrowthrate_=0
			self.areagrowthrate_=0
			return 0
		
		# gr = (log(2*values[-1]+values[-2])-log(2*values[0]+values[1]))/(self.time[-1]-self.time[0])
		# gr = (log(values[-1])-log(values[0]))/(self.time[-1]-self.time[0])
		gr = exp_interpolation((self.time-self.time[0],values),get_model=True)[3]
		if (gr)<=0:
			print("warning: cell {0} has negative {1}growthrate {2:.2f}".format(self.id,attr,gr))
		
		if attr=='length':
			self.lengthgrowthrate_=gr
		elif attr=='maxwidth':
			self.maxwidthgrowthrate_=gr
		elif attr=='area':
			self.areagrowthrate_=gr
		return gr
	
	# does an exponential interpolation on the last num_points timestep.
	# If in the beginning of the cell, takes a portion of length of the parent in the proportion of the length at birth.
	def instantlengthgrowthrateat(self, time, num_points=15):
		celltime = self.time
		celllength=self.length
		try:
			tindex=[np.where(celltime==time)[0][0]]
		except Exception:
			print "Error: cell {} did not exist at time {:.2f}".format(self.id, time)
			return None
		
		if tindex<=num_points-2:
			if self.parent==None:
				return None
			if len(self.parent.time) < (num_points-1)-tindex:
				return None
			celltime=self.parent.time[tindex-(num_points-1):] + celltime[:tindex+1]
			celllength=[l/2 for l in self.parent.length[tindex-(num_points-1):]] + celllength[:tindex+1]
			tindex=(num_points-1)
		
		#TODO: check for possible overflow of a in exp_interpolation if time is too high
		#Is it correctly fixed with -celltime[0]?
		t,exp,a,gr = exp_interpolation((celltime[tindex-(num_points-1):tindex+1]-celltime[0],celllength[tindex-(num_points-1):tindex+1]),get_model=True)
		
		# plt.plot(self.time,self.length,'.k')
		# plt.plot(t,exp)
		
		if (gr)<=0:
			print("warning: cell {} has negative instantaneous lengthsgrowthrate {:.2f}".format(self.id,gr))
			
		return gr

	
	##centered window, window of size num_point
	## truncate window at birth and division
	def __instantlengthgrowthrate(self, num_points=15, recompute=False):

		if hasattr(self, '_instantlengthgrowthrate') and hasattr(self, '_instantlengthgrowthratepoints') and self._instantlengthgrowthratepoints==num_points:
			return self._instantlengthgrowthrate		
		
		
		instantgr = [exp_interpolation((self.time[max(0,t-(num_points)/2):min(len(self.time),t+(num_points+1)/2)],self.length[max(0,t-(num_points)/2):min(len(self.time),t+(num_points+1)/2)]),get_model=True)[3] for t in np.arange(len(self.time))]
		
		self._instantlengthgrowthrate = instantgr
		self._instantlengthgrowthratepoints = num_points
		
		return self._instantlengthgrowthrate
		
	def _parentget(self, attr):
		if self.parent is None:
			return None
		return self.parent.get(attr)
	
	def _sisterget(self, attr):
		if not hasattr(self, 'sibling'):
			if self.parent is None:
				self.sibling = None
			else:
				ss=self.parent.children
				if ss[0].id==self.id:
					self.sibling=ss[1]
					ss[1].sibling=self
				else:
					self.sibling=ss[0]
					ss[0].sibling=self
			
		if self.sibling is None:
			return None
		vv = self.get(attr)
		vs = self.sibling.get(attr)
		if vv==vs==0:
			return None
		return (vv-vs)/(float(vv+vs)/2)
	
	def _daughtersimilarity(self, attr, similarity):
		if len(self.children)!=2:
			return None
		
		if similarity=='dtw':
			return log(dtw_std(moving_average(self.children[0].length),moving_average(self.children[1].length)))
		
		values=[c.get(attr) for c in self.children]
		if None in values:
			return None
		high = np.argmax(values)
		
		poles=[c.get('pole') for c in self.children]
		
		oldpole=None
		if poles[0]==0:
			oldpole=0
		elif poles[1]==0:
			oldpole=1
			
		if similarity=='diff':
			return float(values[high])-values[1-high]
		if similarity=='sum':
			return float(values[high])+values[1-high]
		if similarity=='max':
			return float(values[high])
		if similarity=='min':
			return float(values[1-high])
		if similarity=='avg':
			return (float(values[high])+values[1-high])/2
		elif similarity=='ratio':
			return float(values[high])/values[1-high]
		elif similarity=='oldnewratio':
			if oldpole==None:
				return None
			return float(values[oldpole])/values[1-oldpole]
		elif similarity=='diff':
			return values[high]-values[1-high]
		elif similarity=='oldnewdiff':
			if oldpole==None:
				return None
			return values[oldpole]-values[1-oldpole]
		elif similarity=='percentagediff':
			if np.sum(values) == 0:
				return None
			return float(values[high]-values[1-high])*2/np.sum(values)
		elif similarity=='oldnewpercentagediff':
			if oldpole==None or np.sum(values)==0:
				return None
			return float(values[oldpole]-values[1-oldpole])*2/np.sum(values)
		else:
			print('wrong similarity function: {}'.format(similarity))
			return None
		
	def plot(self):
		plt.figure("cell {}".format(self.get("index")))
		ax=plt.subplot(211, title='length and maxwidth')
		ax.plot(self.time,self.length,'r+')
		ax.plot(self.time,exp_interpolation((self.time-self.time[0],self.length))[1],'r-')
		ax.set_ylabel('length',color='r')
		ax=ax.twinx()
		wmask = self.maxwidth!=None
		ax.plot(self.time[wmask],self.maxwidth[wmask],'b+',)
		ax.plot(self.time[wmask],moving_average(self.maxwidth[wmask],1),'b.-')
		ax.set_ylabel('maxwidth',color='b')
		ax=plt.subplot(212, title='area')
		ax.plot(self.time,self.area,'+')
		ax.plot(self.time,moving_average(self.area,1),'r+-')
		plt.legend()
	
# Lineage class
class Lineage(object):
	
	# docstring TODO
	
	attribute_types={	'id':'nominal',
				'populationgrowthrate':'numeric',
				'areagrowthrate':'numeric',
				'heatshocktime':'numeric',
				'survivalrate':'numeric',
				'populationatshock':'numeric',
				'areaatshock':'numeric',
				'areagrowthrateatshock':'numeric'}
					
	
	def __init__(self, id, experimentid, timestep, fluotimestep=None, cells=[]):
		self.id=id
		self.experimentid=experimentid
		self.cells=cells
		time=np.array([])
		for c in cells:
			time = np.union1d(time,c.time)
		self.time=time
		self.timestep=timestep
		self.fluotimestep=fluotimestep
		
		for c in cells:
			c.lineage=self
		# self._sortlineage()
		
	def __str__(self):
		return "Lineage {}".format(self.id)
		
	def getcells(self):
		return remove_none(self.cells)
	
	def get(self, attribute):
	
		s=attribute.split('_')
		
		if s[0]=='cell':
			return self._cellaggregate('_'.join(s[1:len(s)-1]), s[-1])
		elif attribute=="id":
			return self.id
		elif attribute=="populationgrowthrate":
			if not hasattr(self,"populationgrowthrate"):
				self.populationgrowthrate = self._growthrate("population")
			return self.populationgrowthrate
		elif attribute=="areagrowthrate":
			if not hasattr(self,"areagrowthrate"):
				self.areagrowthrate = self._growthrate("area")
			return self.areagrowthrate
		elif attribute=="heatshocktime":
			if not hasattr(self,"heatshocktime"):
				self.heatshocktime = self._heatshocktime()
			return self.heatshocktime
		elif attribute=="survivalrate":
			if not hasattr(self,"survivalrate"):
				self.survivalrate = self._survivalrate()
			return self.survivalrate
		elif attribute=="populationatshock":
			if self.get("heatshocktime") is None:
				return None
			else:
				return self._populationat(self.get("heatshocktime"))
		elif attribute=="areaatshock":
			if self.get("heatshocktime") is None:
				return None
			else:
				return self._areaat(self.get("heatshocktime"))
		elif attribute=="areagrowthrateatshock":
			if self.get("heatshocktime") is not None:
				return self._areagrowthrateat(self.get("heatshocktime"))
		else:
			print "unvalid attribute for lineage: {}".format(attribute)
			return None
			
	def _sortlineage(self):
		
		self.cells=np.insert(self.cells[np.argsort([c.lefttreeposition for c in remove_none(self.cells)])], 0,None)
		
		i=1
		while (i<len(self.cells)):
			if self.cells[i].lefttreeposition != i:
				self.cells = np.insert(self.cells,i,None)
			i=i+1
		
	def _growthrate(self, type):
	
		if len(remove_none(self.cells)) == 0:
			print "cannot compute growthrate of lineage {} because it contains no cells".format(self.id)
			return None;
		
		if self.get("heatshocktime") is None:
			# the last time (frame) when all cells are recorded
			# we consider the first frame of the last frames of the cells that divide and dont have children
			lastgen = np.array([len(c.children)==0 and c.get("doublingtime") is not None for c in remove_none(self.cells)])
			if not np.any(lastgen):
				print "cannot find when last time is."
				return None
			finaltime = np.min([c.time[-1] for c in remove_none(self.cells)[lastgen]])
		else:#heatshock
			finaltime = self.get("heatshocktime")

		times=np.arange(0,finaltime,self.timestep)
		values=np.zeros(len(times))
		for c in remove_none(self.cells):
			time = np.array(c.time)[np.array([t<finaltime for t in c.time])]
			indexes = [np.where(times==t)[0][0] for t in time]
			for i in indexes:
				v=1 # populationgrowthrate
				if type == "area":
					v=c.areaat(times[i])
				values[i] = values[i] + v
				
		gr = exp_interpolation((times,values),get_model=True)[3]
		if (gr)<=0:
			print("warning: lineage {} has negative {}growthrate ({:.2f})".format(self.id,type,gr))
		
		return gr
		
	def _heatshocktime(self):
		hs=np.array([c.get("alive") is not None for c in self.cells], dtype=bool)
		if not np.any(hs):
			return None
		else:
			# returns the last time of one heatshock cell
			return self.cells[hs][0].time[-1]
	
	
	def _areaat(self, time):
		return np.mean(remove_none([c.areaat(time) for c in self.cells]))
	
	def _areagrowthrateat(self, time):
		
		numframes = 20
		
		if time==0:
			return None
		
		if time/self.timestep < numframes: # to check
			numframes = time/self.timestep
		
		times = np.arange(time-float(numframes)*self.timestep,time+self.timestep, self.timestep)
		validmask = np.zeros(len(self.cells), dtype=bool)
		for i,c in enumerate(self.cells):
			
			if c is None:
				continue
			if c.time[0] > time or c.time[-1] < times[0]:
				continue
			if c.time[-1] < time and len(c.children) < 2 :	#if it divide in the window frams and not all children are here, dont count the cell and remove it's parent if needed
				p=2
				while(validmask[i/p]):	#if the parent was valid
					j=i/p
					q=p
					while(q>0):		# remove descendant starting from the bottom (from i generation)
						for k in range(q*j,q*j+q):
							validmask[k]=False
						q=q/2	# go a generation up when q=1, it is j so stop after
				p=p*2		# check the parent's parent
				continue
			validmask[i] = True
			
		
		# print [c.lefttreeposition for c in self.cells[validmask]]
		
		values = np.empty(len(times))
		for i,t in enumerate(times):
			cellsattmask = np.array([t in c.time for c in self.cells[validmask]])
			
			# print t, cellsattmask
			values[i]=np.sum([c.areaat(t) for c in self.cells[validmask][cellsattmask]])
			
		return exp_interpolation((times,values),get_model=True)[3]
	
	def _populationat(self, time):
		
		return np.sum([time in c.time for c in self.cells])
			
	def _survivalrate(self):
		
		if self.get("heatshocktime") is None:
			return None
			
		values=np.array(remove_none([c.get("alive") for c in self.cells]),dtype=bool)
		return float(np.sum(values))/len(values)
			
	def _cellaggregate(self, attr, aggregate):
	
		values = table(self.cells, [attr])[:,0] # remove_none  ?
		if aggregate in ("avg","mean","average"):
			return np.mean(values)
		else:
			print("unsupported aggregate function: {}".format(aggregate))
			return None

	def __str__(self):
		return "Lineage {}".format(self.id)
			
	def plot(self):
		plt.figure("lineage {}".format(self.id))
		
		ax=plt.subplot(211, title='length')
		# ax.minorticks_on()
		ax.grid(which='both')
		for c in remove_none(self.cells):
			ax.plot(c.time, c.length,'b+' if c.pole==0 else 'r+')
			ax.plot(c.time, moving_average(c.length),'b+-' if c.pole==0 else 'r+-')
			
		ax=plt.subplot(212, title='instantaneous area growthrate', sharex=ax)
		times = np.arange(self.timestep,np.max([np.max(c.time) for c in remove_none(self.cells)]), self.timestep)
		growth = [self._areagrowthrateat(t) for t in times]
		ax.plot(times,growth, 'b+')
		ax.plot(times,moving_average(growth),'b')
		ax.set_ylabel('areagrowthrate',color='b')
		ax=ax.twinx()
		plt.plot(times,[self._populationat(t) for t in times], 'r')
		ax.set_ylabel('#cells',color='r')
		

def getcell(id, cells):

	"""Returns the Cell in cells with the given id"""
	
	cell=cells[np.array([c.id for c in cells])==id]
	if len(cell)>0:
		return cell[0]
	else:
		return None

###################################
# data: select and table                                                 #
################################### 


def select(objects, condition, get_mask=False):

	"""Select a subset array that fits the condition.
	
	Parameters
	-------------
	
	objects: array-like, shape [n_objects,]
		array of Cell or Lineage

	condition: string
		condition is of the form "<attribute> <sign> <value>".
		<attribute> follows the attribute grammar.
		<sign> is "==", "!=", "<", "<=", ">", ">=" or "in".
			if "in", <value> should be a tuple of values.

	get_mask: boolean, default : False

	Returns
	--------

	newarray: numpy.array
		numpy array of Cell or Lineage that fits the condition.
		If get_mask is True, returns mask: numpy array of boolean, shape [n_objects,] such that objects[mask] would normally be returned.

	Examples
	--------

	>>> cells.shape
	(1370,)
	>>> bigcells = select(cells, "avglength > 45")
	>>> bigcells.shape
	(780,)
	>>> wildtype = select(lineages, "id in (746,747,748,749)"
	>>> wildtype[0].id
	746
	>>> mask746_lowgrowth = select(wildtype[0].cells, "lengthgrowthrate <= 0.023", get_mask=True)
	>>> for c in wildtype[0].cells[mask746_lowgrowth]:
	...     print c.get("lengthgrowthrate")
	... 
	0.0173294436084
	0.0202434124752
	0.0215982457257
	0.0228286583693
	
	"""
	
	conditions = condition.split(' && ')  #TODO change to and
	mask = np.ones(len(objects), dtype='bool')
	
	for cond in conditions:
		attribute, symbol, value = _parse_condition(cond, objects[0])
		
		atype = _get_attribute_type(attribute, objects[0])
		
		if attribute == "index": #genealogy condition
			
			if symbol=="descendantof":
				newmask = [descendantof("{}.{}".format(c.lineageid,c.lefttreeposition), value) for c in objects]
			elif symbol=="ancestorof":
				newmask = [ancestorof("{}.{}".format(c.lineageid,c.lefttreeposition), value) for c in objects]
			else:
				print "error parsing genealogy condition: {}".format(cond)
				return
			
		elif(atype == "nominal"): #condition on nominal attribute
			#TODO: check for boolean (=0 or =1)
			
			values = table(objects,[attribute])[:,0]
			
			#TODO : is that efficient? safe?
			newmask = [eval("v {} {}".format(symbol,value if is_float(value) else "'{}'".format(value))) for v in values]
			# newmask = np.empty(len(objects), dtype='bool')
			# for i in xrange(len(objects)):
				# exec("newmask[i] = values[i] %s %s"%(symbol,value))
			
		else: #condition on numerical attribute
			
			values = table(objects,[attribute])[:,0]

			#TODO : is that efficient? safe?
			newmask = [eval("v {} {}".format(symbol,value)) for v in values]
			# newmask = np.empty(len(objects), dtype='bool')
			# for i in xrange(len(objects)):
				# exec("newmask[i] = values[i] %s %s"%(symbol,value))
			
		mask = mask * newmask
	
	if get_mask:
		return mask
	else:
		return objects[mask]

# returns True if value is a descendant of ref(erence)
def descendantof( value, ref):
	v_lid,v_pos = map(int,value.split("."))
	r_lid,r_pos = map(int,ref.split("."))
	if v_lid == r_lid:
		while v_pos>r_pos:
			v_pos = v_pos/2
		return v_pos == r_pos
	return False

# returns True if value is an ancestor of ref(erence)
def ancestorof( value, ref):
	return descendantof( ref, value)
		
	
	
# parse a condition that defines a subset of the cells
# the condition should consist of three parts: the attribute, the symbol, the value(s)
# these should be separated by one (or more) space
# the attribute is an attribute of the cell that can be called with get.
# the symbol can be "==", ">=", ">", "<=", "<", "!=", "in"
# for all except "in" should follow one value.
# for symbol "in", should follow a set of values between brackets
# condition: string
# return a tuple:
#	attribute, symbol, value for numerical/nominal attribute
#	null, symbol, value for genealogy condition (descendant/ancestor)
def _parse_condition(condition, object=None):

	condition = condition.strip()

	if _genealogy_regex.search(condition):
		attribute=None
		symbol, values = condition.strip().split(' ',1)
		values = values.lstrip()
	
	else:
		attribute, s = condition.strip().split(' ',1)
		symbol, values = s.strip().split(' ',1)
		values = values.lstrip()
		
		# checking if condition is valid 
		# compatibility between the different parts
		
		try:
			atype=_get_attribute_type(attribute, object)
		except KeyError:
			raise Exception("Cannot find attribute %s"%(attribute))
		
		number_regex = '(-?[0-9]+(\.[0-9]*)*)'    # does not allow .5 nor scientific notation
		
		_isnumber = re.match(number_regex,values)
		_istuple = re.match('^\(%s(,%s)*\)$'%(number_regex, number_regex),values)
		
		if symbol in ('=','==','!=','>','>=','<','<='):
			if not (_isnumber and atype=="numeric" or atype=="nominal"):
				raise Exception("wrong condition: %s"%(condition))
		elif symbol == 'in':
			if not _istuple:
				raise Exception("wrong condition: %s"%(condition))
		elif symbol in('descendantof','ancestorof'):
			if not atype=="treeindex":
				raise Exception("wrong condition: %s"%(condition))
			
		else:
			raise Exception("Cannot find symbol %s"%(symbol))
			
		
	return attribute,symbol,values
	

def _get_attribute_type(attribute, object=None):
	
	s=attribute.split(' ')

	if len(s)>1: #condition like "attr >= x"
		return 'nominal'

	s=attribute.split('_')

	if s[0]=='parent':
		return _get_attribute_type(attribute[7:], object)
	elif s[0]=='sister':
		return _get_attribute_type(attribute[7:], object)
	elif s[0]=='daughter':
		return _get_attribute_type(s[1], object)
	elif s[0]=='normalised':
		ss=attribute.rsplit('_by_', 1)
		return _get_attribute_type(ss[0][11:], object)
	elif object is not None :
		return object.attribute_types[s[0]]
	else:
		raise Exception("Cannot find attribute of {}".format(attribute))
	
#TODO: modify to use get_attribute_type and do the list of values for nominal attributes
def _get_arff_attribute_type(attribute, object=None):
	
	s=attribute.split(' ')

	if len(s)>1: #condition like "attr >= x"
		return '{0,1}'

	s=attribute.split('_')

	if s[0]=='parent':
		return _get_arff_attribute_type(attribute[7:], object)
	elif s[0]=='sister':
		return _get_attribute_type(attribute[7:], object)
	elif s[0]=='daughter':
		return _get_arff_attribute_type(s[1], object)
	elif s[0]=='normalised':
		ss=attribute.rsplit('_by_', 1)
		return _get_arff_attribute_type(ss[0][11:], object)
	else :
		return _cell_arff_attribute_types[s[0]]
	
	
def table(objects, attributes=[]):
	
	"""Retrieves the attributes of several objects as an array.
	
	Parameters
	-------------
	
	objects: array, shape [n_objects,]
		numpy array of Cell or Lineage

	attributes: list of string
		List of attributes to retrieve the value.
		
		Warning: it has to be a list even for only one attribute.

	Returns
	---------
	
	newarray: numpy.array

		An array where the rows correspond to the objects and the columns are the attributes.

	Examples
	-----------
	
	>>>from subtil.data import table
	>>> table(cells, ["age","lengthgrowthrate"]
	array([[ 0, 0.0273294436084],
			   [ 1, 0.0312488732465],
			   ...,
			   [ 4, 0.0312684213548],
			   [ 3, 0.0264284795423]])
	
	"""
	_data = np.empty((len(objects),len(attributes)), order = 'F', dtype=object)
	
	for ind_a, attr in enumerate(attributes):
				
		s_attr = attr.split(' ')
		if len(s_attr) > 1:  #condition
			mask = select(objects, attr, get_mask=True)
			_data[:, ind_a] = mask.astype('int')
			continue

		s_attr = attr.split('_')
		if s_attr[0] == "normalised":
			ss=attr.rsplit('_by_', 1)
			if len(ss)==1:	# normalise in total
				values = np.append(table(objects, attributes=[attr[11:]]),
									np.ones((len(objects),1), dtype=int),
									axis=1)
			else:	#normalise by some attribute
				values = table(objects, attributes=[ss[0][11:]] + ss[1].split('_and_'))
			
			_uniques = [remove_none(np.unique(v)) for v in values.T[1:]]
			_means = np.empty([len(u) for u in _uniques],dtype=object)
			
			for val in values:
				v = val[0]
				if v is not None:
					# for each cell, get the corresponding tuple of indexes of by_attributes 
					by_attr=tuple([np.where(u==by)[0][0] for u,by in zip(_uniques,val[1:])])
					if _means[by_attr] is None:
						_means[by_attr]=np.array([], float)
					_means[by_attr]=np.append(_means[by_attr],v)
			
			for x in xrange(np.size(_means)):
				if _means.flat[x] is not None:
					_means.flat[x]=_means.flat[x].mean()
					
			for ind_c in range(len(values)):
				v=values[ind_c,0]
				if v is not None:
					by_attr=tuple([np.where(u==by)[0][0] for u,by in zip(_uniques,values[ind_c,1:])])
					_data[ind_c, ind_a] = float(v) / _means[by_attr]
		
		else:
			_data[:, ind_a] = np.array([c.get(attr) for c in objects], order='F', dtype=object)
	
	# _data = np.array([[c.get(attr) for attr in attributes] for c in objects], order='F', dtype=object)
	
	#problem: cannot get parent_normalised attributes TODO
	
	return _data
	
#TO CHECK
def write_arff(cells, attributes, name, file):
	data = np.array(table(cells, attributes),order='C')
	with open(file,'w') as f:
		f.write("@relation %s\n"%(name))
		for a, attr in enumerate(attributes):
			arff_attr = _get_arff_attribute_type(attr)
			if arff_attr[-2:] == '{}':
				arff_attr = "{}{}{}".format(arff_attr[:-1], ','.join(map(str,np.unique(data[:,a]))), "}" )
			f.write("@attribute %s %s\n"%(attr, arff_attr))
		f.write("@data\n")
		for row in data:
			f.write("%s\n"%(','.join(map(str,row)).replace('None','?')))

