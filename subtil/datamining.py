import traceback
import random
import time
from itertools import ifilter
import numpy as np
import pylab as plt
import matplotlib.widgets
import matplotlib.lines
from math import factorial,floor

from sklearn.model_selection import cross_val_predict, GroupKFold
from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import DBSCAN
from sklearn.manifold import MDS

from .utils import remove_none, moving_average
from .data import _get_attribute_type, table, select, descendantof, ancestorof
from .modeltree import ModelTreeRegressor
from .plotting import histplot

#global_variable for subgroups clustering
_sclus=None

# predict the cellcycle at shock of cells which "alive" attribute is not None.
# builds a regression model based on cells that divided
def predictcellcycle(cells, attr=None, regressor=None, normalise=True, evaluate=True, plot_error=False):

	"""Predict cellcycle of heatshocked cells

	Builds a training set to learn a model tree that predicts the cell cycle of heatshocked cells as a number between 0 and 1.

	The filamentous cells, defined as cells that get bigger than 100px, are ignored in the construction of the model.

	Parameters
	----------
	
	cells: array of cells
		This array contains both heatshock cells and fully divided cells from which the regression model is learned.
		
	attr: list of string
		List containing the attributes to use for the prediction.
		Default list is ["generation", "age", "youth", "parent_doublingtime", "parent_lengthincrease", "parent_lengthgrowthrate",  "lengthatbirth", "lengthgrowthrate","lengthgrowthrateatshock", "timeatshock", "lengthincrease", "absolutefluo", "absfluoincrease"]

	regressor: sklearn.Regressor
		The regressor to use. Default is the ModelTreeRegressor.

	normalise: boolean, default: True
		Whether to normalise the attributes. If True, each attribute is normalised between 0 and 1 before learning the model.
	
	evaluate: boolean, default: True
		Whether to perform 10-fold cross validation.

	plot_error: boolean, default: False
		Whether to plot the prediction error. Ignored if evaluate is False.

	Returns
	-------
	

	"""

	hs = np.array([c.get("alive") is not None for c in cells]) #heatshock cells
	
	if attr is None:
		attr=["generation","age","youth","parent_doublingtime","parent_lengthincrease","parent_lengthgrowthrate", "lengthatbirth","lengthgrowthrate","lengthgrowthrateatshock","timeatshock", "lengthincrease","absolutefluo", "absfluoincrease"]
	
	# building training set.
	# excluding:
	#	cells that don't divide
	#	cells of first generation
	#	cells whose fluo has not been measured if fluo attributes are required
	train_X = []
	train_y = []
	linid = []
	for c in cells[hs==False]:
		#we need cell that divide
		if c.get("doublingtime") is None:
			continue
		if c.get("generation")<1:
			continue
		if (c.lineage.fluotimestep is None or c.lineage.fluotimestep == 0) and np.any(["absolutefluo" in attr,"absfluoincrease" in attr]):
			raise Exception("cell {} has no fluorescence")
		# choosing t so that is on a frame where fluo was actually measured (every 8 frames, i.e. every 4s)
		indexes = range(len(c.time))
		if (c.lineage.fluotimestep is not None and c.lineage.fluotimestep != 0) and np.any(["absolutefluo" in attr,"absfluoincrease" in attr]):
			indexes=np.where((np.array(c.time)%c.lineage.fluotimestep)==0)[0]
		if len(indexes)==0:
			continue
		t= random.choice(indexes)
		time=c.time[t]-c.get("birthtime")
		lengthatbirth=c.get("lengthatbirth")
		length=c.length[t]
		# length=c.length[t]
		if length>100: #ignoring cells that became filamentous
			continue
		# fluoatbirth = c.fluorescence[0]
		# fluoatbirth = c.fluorescence[0]*c.area[0]
		# fluo = c.fluorescence[t]
		if "absolutefluo" in attr:
			fluo = c.totalfluorescenceat(c.time[t])
		if "absfluoincrease" in attr:
			fluoincrease = c.totalfluorescenceat(c.time[t])-c.totalfluorescenceat(c.time[indexes[0]])
		inst_gr = c.instantlengthgrowthrateat(c.time[t])
		if "lengthgrowthrateatshock" in attr and inst_gr == None:
			continue
			# ignoring inst_gr because it is None for some cells
		X=[]
		for a in ["generation","age","youth","parent_doublingtime","parent_lengthincrease","parent_lengthgrowthrate", "lengthatbirth", "lengthgrowthrate"]:
			if a in attr:
				X.append(c.get(a))
		if "timeatshock" in attr:
			X.append(time)
		if "lengthincrease" in attr:
			X.append(length-lengthatbirth)
		if "absolutefluo" in attr:
			X.append(fluo)
		if "absfluoincrease" in attr:
			X.append(fluoincrease)
		if "lengthgrowthrateatshock" in attr:
			X.append(inst_gr)
			
		train_X.append(X)
		# train_X.append([time, length-lengthatbirth, fluoincrease])
		train_y.append(float(time)/c.get("doublingtime"))

		linid.append(c.get("lineageid"))
			
	if normalise:
		scaler = MinMaxScaler()
		train_X = scaler.fit_transform(np.array(train_X)).tolist()
	else:
		scaler = None
	
	# return train_X,train_y
	
	if regressor is None:
		regressor = ModelTreeRegressor()
	
	if evaluate:
		predicted_y=0
		try:
			#TODO: change LabelKFold to GroupKFold 
			lkf = LabelKFold(linid, n_folds=10)
			
			predicted_y=cross_val_predict(regressor, np.array(train_X), np.array(train_y), cv=lkf)
			regr_mse=((predicted_y-np.array(train_y))**2).sum()
			res_mse=((np.array(train_y)-np.mean(train_y))**2).sum()
			perf=1-regr_mse/res_mse
			if plot_error:
				plt.figure(str(regressor))
				plt.scatter(np.array(train_y),predicted_y)
				plt.xlabel("cycle")
				plt.ylabel("predictedcycle")
				plt.xlim([0,1])
				plt.ylim([0,1])
		except Exception:
			print("Couldn't evaluate given regressor.")
			print(traceback.format_exc())
			perf=0
	else:
		perf=0
	
	
	# regressor.fit(train_X,train_y)
	regressor.fit(np.array(train_X),np.array(train_y))
	
	for c in cells[hs]:
		try:
			t=len(c.time)-1
			time=c.time[t]-c.get("birthtime")
			lengthatbirth=c.get("lengthatbirth")
			length=c.length[t]
			if length>100: #ignoring filamentous cells
				continue
			if "absolutefluo" in attr:
				fluo = c.totalfluorescenceat(c.time[t])
			if "absfluoincrease" in attr:
				fluoincrease = c.totalfluorescenceat(c.time[t])-c.totalfluorescenceat(c.lineage.fluotimestep*floor(1+(c.time[0]-c.lineage.timestep)/c.lineage.fluotimestep) )
			inst_gr = c.instantlengthgrowthrateat(c.time[t])
			predict_X=[]
			for a in ["generation","age","youth","parent_doublingtime","parent_lengthincrease","parent_lengthgrowthrate", "lengthatbirth","lengthgrowthrate"]:
				if a in attr:
					predict_X.append(c.get(a))
			if "timeatshock" in attr:
				predict_X.append(time)
			if "lengthincrease" in attr:
				predict_X.append(length-lengthatbirth)
			if "absolutefluo" in attr:
				predict_X.append(fluo)
			if "absfluoincrease" in attr:
				predict_X.append(fluoincrease)
			if "lengthgrowthrateatshock" in attr:
				predict_X.append(inst_gr)
			# predict_X=[[c.generation, c.age, c.youth, c.get("parent_doublingtime"), c.get("parent_lengthincrease"), c.get("parent_lengthgrowthrate"), lengthatbirth,time,  length-lengthatbirth, fluo, fluoincrease]]
			# predict_X=[time, length-lengthatbirth, fluoincrease]
			
	
			if normalise:
				predict_X = scaler.transform(np.array([predict_X])).tolist()

			predictedcycle = regressor.predict(predict_X)[0]
			if predictedcycle<0:
				predictedcycke=0
			if predictedcycle>1:
				predictedcycle=1
			c.cycleatshock = predictedcycle

		except Exception:
			print("error when computing cycleatshock of cell {}".format(c.id))
			print(traceback.format_exc())
			c.cycleatshock=None
			
	Cell.attribute_types['cycleatshock']='numeric'

	if evaluate:
		return perf, train_X, train_y, predicted_y, scaler
	else:
		return perf, train_X, train_y, scaler

	
	
def search(objects, description_attributes, target, depth=2, k=100, minsize=20, beamsize=2000):
	
	"""Performs subgroup discovery.

	Looks for subgroups defined by the attributes in description_attributes with the defined target.

	Uses beam search to look for subgroups.
	
	The quality measure is picked automatically depending on the type of the target value.
	

	Parameters
	----------
	
	objects: array-like
		This array contains objects (Cell, Lineage) to be mined.
		Subgroups will consists of subset of this array.
		
	description_attributes: list of string
		List containing the attributes to use as description attributes.

	target: string
		The target attribute, possibly preceded by "high" or "low" for numerical attributes
	
	depth: int, default: 2
		The depth of the beam search.

	k: int, default: 100
		Number of subgroups to return.

	minsize: int, default:0
		Minimum number of instances.
		If set to 0 (default), the minimal size is equal to 1 percent of the total number of objects

	beamsize: int, default:2000
		Size of the beam to use in the beam search.

	Returns
	-------
	subgroups: numpy array
		An array of Subgroup objects.
	
	"""



	assert depth > 0, "Depth cannot be less than 1."

	if minsize==0:
		minsize = ceil(len(objects/100))
		
	beamsize=beamsize
		
	beam = np.array([Subgroup(objects,
					mask = np.ones(len(objects),dtype = bool), 
					description = "",
					target = target,
					quality = 0)])
	subgroups = beam
	
	for d in range(1,depth+1):
		candidates = beam
		progress = 0
		print("depth {}".format(d))
		start=time.time()
		print("refining...")
		for i,sg in enumerate(beam):
			for attr in description_attributes:
				try:
					candidates = np.append(candidates, _filter_size(_refine(sg, attr),minsize))
				except Exception as ex:
					print "Error while refining with attribute %s"%(attr)
					print(traceback.format_exc())
				
		end=time.time()
		print("    {} s".format(end-start))
		
		print("computing quality...")
		start=time.time()
		_compute_quality_subgroups(candidates, target=target)
		end=time.time()
		print("    {} s".format(end-start))
		
		print("filtering...")
		start=time.time()
		beam=_filter_cover(candidates, th=0.5, beamsize=beamsize)
		subgroups = _filter_cover_double(subgroups,beam, th=0.5, beamsize=k)
		end=time.time()
		print("    {} s".format(end-start))
		
	
	subgroups = _reduce_descr(subgroups)

	print("done.")
	
	return subgroups

# return subgroups of that subgroup
def _refine(subgroup, attribute):

	#print("refining {0} with {1}...".format(subgroup.description, attribute))
	
	atype = _get_attribute_type(attribute, subgroup.cells[0])
	
	values = table(subgroup.cells[subgroup.mask], [attribute])[:,0]
	
	subgroups=[subgroup]
	
	if atype == "numeric":
		numbins = 8
		u_values = remove_none(np.unique(values))
		if len(u_values) <= numbins:
			bins = u_values
		else:
			bins = np.unique(np.percentile(remove_none(values),np.append(np.arange(0,100,100/numbins),100).tolist()))
		
		
		for bin in bins[1:-1]:
			infmask = np.array(subgroup.mask)
			supmask = np.array(subgroup.mask)
			#j=0
			#for i,m in enumerate(subgroup.mask):
			#	if m:
			#		if values[j]<bin:
			#			supmask[i]=False
			#		else:
			#			infmask[i]=False
			#		j=j+1
			infmask[subgroup.mask] = infmask[subgroup.mask]*(values<bin)
			supmask[subgroup.mask] = supmask[subgroup.mask]*(values>=bin)
			
			if not np.all(subgroup.mask==infmask):
				subgroups.append(Subgroup(cells = subgroup.cells,
					mask = infmask, 
					description = subgroup.description + (" && " if subgroup.description != "" else "") + "{0} < {1}".format(attribute,bin),
					target = subgroup.target,
					quality = None))
			
			if not np.all(subgroup.mask==supmask):
				subgroups.append(Subgroup(cells = subgroup.cells,
					mask = supmask, 
					description = subgroup.description + (" && " if subgroup.description != "" else "") + "{0} >= {1}".format(attribute,bin),
					target = subgroup.target,
					quality = None))
					
			
	elif atype == "nominal":
		u_values = remove_none(np.unique(values))
		for v in u_values:
			mask=np.array(subgroup.mask)
			#j=0
			#for i,m in enumerate(mask):
			#	if m:
			#		if values[j]!=v:
			#			mask[i]=False
			#		j=j+1
			mask[subgroup.mask] = mask[subgroup.mask]*(values==v)
					
			subgroups.append(Subgroup(cells = subgroup.cells,
					mask = mask, 
					description = subgroup.description + (" && " if subgroup.description != "" else "") + "{0} == {1}".format(attribute,v),
					target = subgroup.target,
					quality = None))
			
		# descriptions = ["%s == %f"%(description_attribute,v) for v in u_values] + ["%s != %f"%(description_attribute,v) for v in u_values]
	elif atype == "list":
		pass
	elif atype == "treeindex":
		descdone = []
		for v in values:
			#ancestorof and descendantof v
			ancmask=np.array(subgroup.mask)
			descmask=np.array(subgroup.mask)
			j=0
			for i,m in enumerate(ancmask):
				if m:
					if not ancestorof(values[j],v):
						ancmask[i]=False
					if not descendantof(values[j],v):   #elif?
						descmask[i]=False
					j=j+1
			
			descdone.append(v)
					
			if not np.all(subgroup.mask==ancmask):
				subgroups.append(Subgroup(cells = subgroup.cells,
					mask = ancmask, 
					description = subgroup.description + (" && " if subgroup.description != "" else "") + "{0} ancestorof {1}".format(attribute,v),
					target = subgroup.target,
					quality = None))
			if not np.all(subgroup.mask==descmask):
				subgroups.append(Subgroup(cells = subgroup.cells,
					mask = descmask, 
					description = subgroup.description + (" && " if subgroup.description != "" else "") + "{0} descendantof {1}".format(attribute,v),
					target = subgroup.target,
					quality = None))
					
			#descendantof ancestors of v
			for a in _ancestor_indexes(v):
				if a in descdone:
					continue
				mask=np.array(subgroup.mask)
				j=0
				for i,m in enumerate(mask):
					if m:
						if not descendantof(values[j],a):
							mask[i]=False
						j=j+1
						
				if not np.all(subgroup.mask==mask):
					subgroups.append(Subgroup(cells = subgroup.cells,
						mask = mask, 
						description = subgroup.description + (" && " if subgroup.description != "" else "") + "{0} descendantof {1}".format(attribute,a),
						target = subgroup.target,
						quality = None))
				descdone.append(a)
	
	return np.array(subgroups)
	
def _ancestor_indexes(index):
	lin,pos=map(int,index.split('.'))
	anc=[]
	while pos>1:
		pos=pos/2
		anc=anc+[pos]
	return ["{0}.{1}".format(lin,p) for p in anc]
	
###########
## removes subgroups of size less than minsize
###############
def _filter_size(subgroups, minsize=10):
	if len(subgroups)==0:
		return subgroups
	return subgroups[np.array([s.size>=minsize for s in subgroups])]
	
###################
# remove unnecessary conditions from subgroups description
# ex: if description= "A>5 && A>7", it changes it to "A>7"
###################
def _reduce_descr(subgroups):
	if len(subgroups)==0:
		return subgroups
	for subgroup in subgroups:
	
		# print ""
		# print subgroup.description
	
		predicates=subgroup.description.split(" && ")
		i=0
		while(i<len(predicates)-1):
			j=i+1
			while(j<len(predicates)):
				if predicates[i]==predicates[j]:
					predicates=predicates[:i]+predicates[i+1:]
					i=i-1
					break
				descr=predicates[:i]+predicates[i+1:j]+predicates[j+1:]
				cells=subgroup.cells if len(descr)==0 else select(subgroup.cells, " && ".join(descr))
				maski=select(cells,predicates[i],get_mask=True)
				maskj=select(cells,predicates[j],get_mask=True)
				if np.sum(maski[maskj==False])==0:
					predicates=predicates[:j]+predicates[j+1:]
				elif np.sum(maskj[maski==False])==0:
					predicates=predicates[:i]+predicates[i+1:]
					i=i-1
					break
				else:
					j=j+1
			i=i+1
		newdescr=" && ".join(predicates)
		# if newdescr != subgroup.description:
			# print "{} -> {}".format(subgroup.description, newdescr)
		subgroup.description = newdescr
		# print subgroup.description
	return subgroups
	
###########################
# remove redundant subgroup (i.e. subgroups with similar cover)
# deprecated (correct but slower)
##########################
def _filter_cover_old(subgroups, th=0.33, beamsize = 1000):
	if len(subgroups)==0:
		return subgroups
	
	subgroups = (subgroups[np.argsort([s.quality for s in subgroups])])[::-1].tolist()
	
	for i,sg1 in enumerate(subgroups):
		if i>beamsize:
			break
		if np.any([_cover_dissimilarity(sg1,sg2) < th for sg2 in subgroups[:i]]):
			subgroups.remove(sg1)

	return np.array(subgroups[:i])

###########################
# remove redundant subgroup (i.e. subgroups with similar cover)
##########################
def _filter_cover(subgroups, th=0.33, beamsize = 1000):
	if len(subgroups)==0:
		return subgroups

	subgroups = (subgroups[np.argsort([s.quality for s in subgroups])])[::-1].tolist()
	
	i=0
	for i,sg1 in enumerate(subgroups):
		if i == beamsize:
			break
		for sg2 in ifilter(lambda sg2: _cover_dissimilarity(sg1,sg2) < th, subgroups[i+1:]):
			subgroups.remove(sg2)

	return np.array(subgroups[:i+1])

def _filter_cover_double(subgroups1, subgroups2, th=0.33, beamsize = 1000):
	if len(subgroups1)==0 and len(subgroups2)==0:
		return subgroups1
	mask1=np.zeros(len(subgroups1),dtype="bool")
	mask2=np.zeros(len(subgroups2),dtype="bool")
	i=0
	j=0
	n=0
	m=0
	
	#subgroups1 = (subgroups1[np.argsort([s.quality for s in subgroups1])])[::-1]
	#subgroups2 = (subgroups2[np.argsort([s.quality for s in subgroups2])])[::-1]
	
	while(n+m<beamsize and i<len(subgroups1) and j<len(subgroups2)):
		sg1 = subgroups1[i]
		sg2 = subgroups2[j]
		if sg1.quality >= sg2.quality:
			k=0
			while(k<m):
				sg3 = subgroups2[k]
				diss = _cover_dissimilarity(sg1,sg3)
				if diss<th:
					break
				k=k+1
			if k==m:
				mask1[i]=True
				n=n+1
			i=i+1
		else:
			k=0
			while(k<n):
				sg3 = subgroups1[k]
				diss = _cover_dissimilarity(sg2,sg3)
				if diss<th:
					break
				k=k+1
			if k==n:
				mask2[j]=True
				m=m+1
			j=j+1

	subgroups = np.append(subgroups1[mask1],subgroups2[mask2])
	subgroups = (subgroups[np.argsort([s.quality for s in subgroups])])[::-1]

	return subgroups

def _compute_quality_subgroups(subgroups, target=""):
	
	if len(subgroups) == 0:
		return
	
	if target=="":
		target = subgroups[0].target

	if "low" in target or "high" in target:
		target_direction = target.split(" ")[0]
		target_attribute = " ".join(target.split(" ")[1:])
	else:
		target_attribute = target
		target_direction = "any"

	target_type = _get_attribute_type(target_attribute, subgroups[0].cells[0])

	if target_type == "numeric" or target_type == "nominal":
	
		values = table(subgroups[0].cells, [target_attribute])[:,0]
	
		for subgroup in subgroups:
			subgroup.quality = _compute_quality(values, subgroup.mask, target_type, target_direction)

	elif targettype == "list":
		
		avg = average_cell(subgroups[0].cells)

def _compute_quality(values, mask, target_type, target_direction="any"):
	
	if target_type == "numeric":
		#uses WRAcc: Weighted Z-score
		
		muD = np.mean(remove_none(values))	#average dataset
		sigmaD = np.std(remove_none(values))	#std deviation dataset
		D = len(remove_none(values))			 			#cardinality dataset
			
		muS = np.mean(remove_none(values[mask]))		#average subgroup
		# S = np.sum(subgroup.mask)					#cardinality subgroup
		S = len(remove_none(values[mask]))			#cardinality subgroup
			
			
#		measure = (float(muS)/muD - 1)*pow(float(S)/D,0.5)*100
		measure = ((float(muS)-muD)/sigmaD)*pow(float(S)/D,0.5)*100
			
		# print muD, muS, D, S, measure
		
		if target_direction == "high":
			return measure
		elif target_direction == "low":
			return -measure
		elif target_direction == "any":
			return abs(measure)
				
	elif target_type == "nominal":
		
		u_values = np.unique(remove_none(values))

		if len(u_values)==2 and True in u_values and False in u_values: #Boolean : use Precision
			nD=np.sum(values)
		
			# mD, nD = mode(remove_none(values))	#mode and count
			D = len(remove_none(values))			#cardinality dataset
		
			nS=np.sum(values[mask])
			# nS = np.sum(remove_none(values[subgroup.mask])==mD)
			S = len(remove_none(values[mask]))		#cardinality subgroup
	
			if S==0:
				measure=0
			else:
				# measure = (nS/S-nD/D)*pow(float(S)/D,0.5)*100
				# measure = np.sum([abs(float(s)/S-float(d)/D) for s,d in nS,nD])*pow(float(S)/D,1)
				measure = (float(nS)/S - float(nD)/D)*pow(float(S)/D,0.5)
				# measure = np.sum([d-s for d,s in nD,nS])*pow(float(S)/D,0.5)
#			print nD, nS, D, S, measure
			
			return measure

		else: #Multi-class
		#uses MWRAcc: Multi-class Weighted Relative Accuracy
		
			nD=np.array([np.sum(values==v) for v in u_values])
		
			# mD, nD = mode(remove_none(values))	#mode and count
			D = len(remove_none(values))			#cardinality dataset
		
			nS=np.array([np.sum(values[mask]==v) for v in u_values])
			# nS = np.sum(remove_none(values[subgroup.mask])==mD)
			S = len(remove_none(values[mask]))		#cardinality subgroup
	
			if S==0:
				measure=0
			else:
				# measure = (nS/S-nD/D)*pow(float(S)/D,0.5)*100
				# measure = np.sum([abs(float(s)/S-float(d)/D) for s,d in nS,nD])*pow(float(S)/D,1)
				measure = (abs(nS.astype(float)/S - nD.astype(float)/D).sum())*pow(float(S)/D,0.5)
				# measure = np.sum([d-s for d,s in nD,nS])*pow(float(S)/D,0.5)
	
			#if target == "high":
			#	return measure
			#elif target == "low":
			#	return -measure
			#elif target == "any":
			#	return abs(measure)
			return measure
				
	elif target_type == "list":
		pass
		#TODO
	
####################
# keep only significant subgroups
# draws n random subgroup, compute their quality
# keep subgroups with experimental pvalue at least p
def _filter_significant(subgroups, p=0.01, n=200):
	
	if len(subgroups) == 0:
		return subgroups
	
	if " " in subgroups[0].target:
		target_attribute, target_direction = subgroups[0].target.split(" ")
	else:
		target_attribute = subgroups[0].target
		target_direction = "any"

	target_type = _get_attribute_type(target_attribute, subgroups[0].cells[0])

	#maxsize=len(subgroups[0].cells)
	values = table(subgroups[0].cells, [target_attribute])[:,0]

	if target_type=="nominal" and len(np.unique(values))==2: #Binary

		pbinom = np.mean(values)
		sizes = np.array([sg.size for sg in subgroups])
		ks = np.array([values[sg.mask].sum() for sg in subgroups])
		
		probas = np.array([(factorial(size)/(factorial(k)*factorial(size-k)))*(pbinom**(k))*((1-pbinom)**(size-k)) for size,k in zip(sizes, ks)]) #binomial (should be hypergeometric, but we assume len(cell) >> sg.size
		return subgroups[probas<p]

	qualities = np.empty((len(subgroups),n))
	for i in xrange(n):
		start = time.time()
		random.shuffle(values)
		#nonemask = np.array([v is not None for v in values])
		 #print "shuffling:{}".format(time.time()-start)
		#qualities[:,i] = [_compute_quality(values, sg.mask, target_type, "any") for sg in subgroups]
		qualities[:,i] = [_compute_quality(values, sg.mask, target_type, target_direction) for sg in subgroups]
		#print "compute_qualities:{}".format(time.time()-start)
		#print time.time()-start
	thresholds = np.percentile(qualities, (1-p)*100, axis=1)
	
	return subgroups[np.array([sg.quality>th for sg,th in zip(subgroups, thresholds)])]
	
class Subgroup(object):

	"""TODO"""

	def __init__(self, cells, mask, description, target, quality = 0):
		self.cells = cells
		self.mask = mask
		self.description = description
		self.target = target
		self.quality = quality
		self.size=np.sum(mask)
			
	def __str__(self): return "Subgroup, description:{}, target:{}, quality:{:0.2f}".format(self.description,self.target,self.quality if self.quality is not None else float('NaN'))
	
	def plot(self):
		plt.figure(self.description)
		target_attribute = self.target.split(" ")[-1]
		#values=table(self.cells,[target_attribute])[:,0]
		#x,bins,patches=plt.hist(remove_none(values), bins = 15)
		#plt.hist(remove_none(values[self.mask]),bins=bins)
		#plt.xlabel(target_attribute)
		histplot(self.cells, target_attribute, yseries=self.mask, width=1)


def _cover_dissimilarity(sg1, sg2):
		diss=float(np.sum(sg1.mask^sg2.mask))/max(np.sum(sg1.mask),np.sum(sg2.mask))
		return diss



def cluster_subgroups(subgroups):

	"""TODO"""

	global _sclus
	_sclus = SubgroupClustering()
	_sclus.fit(subgroups)	
	

class SubgroupClustering(object):
	
	def __init__(self):
		self.eps = 0.57
		self.ax = None
		self.subgroups = None
		self.distances = None
		self.artist = None
		self.artists = []

	def fit(self, subgroups):

		self.subgroups = subgroups
		self.distances = np.array([[_cover_dissimilarity(s1, s2) for s2 in subgroups] for s1 in subgroups])
		
		self.eps=np.percentile(self.distances, (100+3*(len(subgroups)-1))/len(subgroups))
		
		clus = DBSCAN(eps=self.eps, min_samples=1, metric="precomputed")
		clus.fit(self.distances)

		manifold = MDS(dissimilarity="precomputed")
		#manif = TSNE(metric="precomputed")

		self.coord = manifold.fit_transform(self.distances)
		cm = plt.get_cmap()
		colors = np.array([cm(i) for i in np.arange(0,1,1./(max(clus.labels_)+1))])
		np.random.shuffle(colors)

		maxq = np.max([sg.quality for sg in subgroups]) + 0.0001
		minq = np.min([sg.quality for sg in subgroups]) + 0.0001
	
		fig= plt.figure()
		self.ax = plt.axes([0.1,0.2,0.8,0.7])

		plt.subplots_adjust(bottom=0.2)
	
		for i,sg in enumerate(self.subgroups):
			plt.plot(self.coord[i,0],self.coord[i,1], c = colors[clus.labels_[i]], label = sg.description, marker = 'o', markersize= 5+25*(sg.quality-minq)/(maxq-minq), markeredgecolor=(0,0,0,0), picker=5)

		plt.xticks(())
		plt.yticks(())

		fig.canvas.mpl_connect('button_press_event', self.click)
		fig.canvas.mpl_connect('motion_notify_event', self.hover)

		self.epsax = plt.axes([0.2,0.1,0.6,0.03])
		self.slider = matplotlib.widgets.Slider(self.epsax, 'eps', valmin=0.0, valmax=1.0, valinit=self.eps, valfmt='%0.2f')
		self.slider.on_changed(self.slide)

	def slide(self, eps):

		self.eps=eps
	
		clus = DBSCAN(eps=eps, min_samples=1, metric="precomputed")
		clus.fit(self.distances)

		cm = plt.get_cmap()
		colors = np.array([cm(i) for i in np.arange(0,1,1./(max(clus.labels_)+1))])
		np.random.shuffle(colors)

		for line,c in zip(self.ax.findobj(matplotlib.lines.Line2D),clus.labels_):
			line.set_c(colors[c])
		plt.draw()

	def click(self, event):

		if event.inaxes == self.ax:

			for art in self.artists[::-1]:
				if art.contains(event)[0]:
					art.remove()
					plt.draw()
					self.artists.remove(art)
					return

			if event.dblclick:
				for sg in self.subgroups:
					if self.artist is not None and self.artist.get_text() == sg.description:
						sg.plot()
						return

			if self.artist is not None:
				for art in self.artists:
					if art.get_text() == self.artist.get_text():
						return
				self.artists.append(self.ax.text(self.artist.get_position()[0], self.artist.get_position()[1], self.artist.get_text(), bbox=dict(facecolor='white', alpha=1.)))

	def hover(self, event):

		if event.inaxes == self.ax:
		
			if self.artist is not None:
				self.artist.remove()
				self.artist = None

			for art in self.artists:
				if art.contains(event)[0]:
					plt.draw()
					return

			ind=None
			ex = event.xdata
			ey = event.ydata
	
			lines = self.ax.findobj(matplotlib.lines.Line2D)
			dist = np.inf
			ind = -1
			for i,line in enumerate(lines):
				if line.contains(event):
					newdist = ((line.get_xdata() - event.xdata)**2 + (line.get_ydata() - event.ydata)**2)**0.5
					if newdist<dist:
						ind = i
						dist = newdist
			if ind>=0:
				self.artist = self.ax.text(lines[ind].get_xdata(), lines[ind].get_ydata(), lines[ind].get_label(), bbox=dict(facecolor='white', alpha=1.))
				plt.draw()
			else:
				plt.draw()
	
############################# Extra #############################
# Stuff for time-series sd. remove or leave it?
def average_cell(cells):

	"""Computes an average length time-series

	Warning: all cells must have same timestep, and no missing values (That last point should be fixed).

	Parameters
	-------------
	
	cells: array of cells
		The cells to average.
		
	Returns
	---------

	average_cell: array-like
		The length series of the average cell.

	"""

	#init: take the longest cells and shifting is 0.
	refcell = moving_average(cells[np.argmax([len(c.length) for c in cells])].length)
	shifts = np.zeros(len(cells))
	cost = np.sum([_align_cost(refcell,c.length,0) for c in cells])

	it=0
	while True:

		#shifting
		new_shifts = np.array([_align_seq(refcell, c.length) for c in cells])
		new_shifts = new_shifts-int(round(np.mean(new_shifts)))

		#averaging
		new_refcell = _ref_cell(cells, new_shifts)

		new_cost = np.sum([_align_cost(new_refcell,c.length,s) for c,s in zip(cells, new_shifts)])

		if new_cost >= cost or it >= 20:
			return refcell[:-1]
		else:
			refcell = new_refcell
			shifts = new_shifts
			cost= new_cost
			it = it+1
	
def _align_seq(sequence1, sequence2):

	#TODO: improve!!! Way too slow!

	shift=np.argmin([_align_cost(sequence1,sequence2,s) for s in range(-len(sequence2)+1,len(sequence1))])-len(sequence2)+1

	return shift
	
def _align_cost(seq1, seq2, shift):

	cost=0
	i=0
	for x,y in (zip(seq1[shift:],seq2) if shift>=0 else zip(seq1,seq2[-shift:])):
		cost=cost+pow(float(x-y),2)
		i=i+1

	if i==0:
		if shift < 0:
			return pow(float(seq2[-1]-seq1[0]),2)
		else:
			return pow(float(seq2[0]-seq1[-1]),2)
	else:
		return cost/i


def _ref_cell(cells, shifts):

	#mean of shifts should be 0.

	reflen = int(round(np.mean([len(c.time)+s for c,s in zip(cells,shifts)])))
	avg = np.zeros(reflen, dtype=float)
	count = np.zeros(reflen, dtype=float)

	for c,s in zip(cells,shifts):
		if s>reflen or s<-len(c.time):
			continue
		avg[max(0,s):min(reflen,len(c.time)+s)] = avg[max(0,s):min(reflen,len(c.time)+s)] + c.length[max(0,-s):min(len(c.time),reflen-s)]
		count[max(0,s):min(reflen,len(c.time)+s)] = count[max(0,s):min(reflen,len(c.time)+s)] + 1
		
	return avg/count

def _plotavg(cells, shifts, refcell):

	#mean(shifts) should be 0.

	timestep = cells[0].time[1]-cells[0].time[0]
	for c,s in zip(cells,shifts):
		plt.plot(c.time-c.time[0]+s*timestep, c.length, color=(0.5,0.5,0.5,0.5))
	plt.plot(np.arange(len(refcell))*timestep, refcell, linewidth=3)
	plt.draw()





