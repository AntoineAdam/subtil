import traceback
import numpy as np
import pylab as pl 
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
from math import log,sqrt,exp,isnan
from sklearn import linear_model

from . import data

	
#position : "lefttreeposition" or "outtreeposition"
def treeplot(lineages, size, color, size_missing = "avg", color_range="auto", color_min = (1,0,0), color_max=(0,0,1), color_missing=(0.5,0.5,0.5), linewidth=5, position="lefttreeposition", average=False):

	"""TODO"""
	
	if average:
		#TODO
		print "average not supported yet"
		pass
	else:
		for lin in lineages:
			cells = lin.cells
			mask = data.select(cells, "generation <= 12", get_mask=True)
			if np.any(mask==False):
				print "Warning: The tree is too big to be displayed and was truncated beyond generation 12"
				cells = cells[mask]
			vdata = data.table(cells, [position, size, color])
			vdata = vdata[np.argsort(vdata[:,0],axis=0)]

			values = np.zeros((np.nanmax(vdata[:,0])+1,2))*np.nan
			for p,s,c in vdata:
				values[p,:] = [s,c]

			if color_range=="auto":
				_color_range=[np.nanmin(values[:,1]),np.nanmax(values[:,1])]
			else:
				_colorc_range=color_range
			if size_missing == "avg":
				_size_missing=np.nanmean(values[:,0])
			else:
				_size_missing=size_missing

			#print values

			plt.figure("lineage {}".format(lin.id))

			ax = plt.gca()
			cax,cbargs = mpl.colorbar.make_axes(ax, location='right', pad=0.05, fraction = 0.07)
			
			plt.sca(ax)
			_drawtree(values[1:], _color_range, color_min, color_max, color_missing, _size_missing, linewidth)
			plt.ylabel(size)
			
			norm = mpl.colors.Normalize(vmin = _color_range[0], vmax = _color_range[1])
			cmap = mpl.colors.LinearSegmentedColormap.from_list("custommap", [color_min, color_max], N=len(np.unique(values[:,1])))
			mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation = cbargs["orientation"], ticklocation = cbargs["ticklocation"], label = color)

			
	
def _drawtree(values, color_range, color_min, color_max, color_missing, size_missing, linewidth):
	#values is a nx2 table where columns are (size, color)
	#values has to be of length 1, 3, 7, 15, ... 2^n-1?????????

	gen = int(log(len(values), 2))

	ax = plt.gca()
	plt.xticks([])
	#plt.xlim([float(-1*vscale)/scale,float(pow(2,gen)*vscale)/scale])
	#plt.ylim([float((-gen-2)*hscale)/scale, 0])

	#computing x coordinates
	coords = np.zeros((2**(gen+1),2))
	coords[2**gen:,0] = np.arange(pow(2,gen), dtype=float)
	for i in range(2**gen-1,0,-1):
		coords[i,0]=np.mean(coords[2*i:2*(i+1),0])
	coords = coords[1:]

	#computing y coordinates
	for i in range(1, len(values)):
		if np.isnan(values[(i+1)/2-1,0]):
			coords[i,1]=coords[(i+1)/2-1,1]+size_missing
		else:
			coords[i,1]=coords[(i+1)/2-1,1]+values[(i+1)/2-1,0]
		
		ax.add_line(mpl.lines.Line2D([coords[i,0],coords[(i+1)/2-1,0]], [coords[i,1],coords[i,1]], color=(0.5,0.5,0.5), zorder=1))

	minv = color_range[0]
	maxv = color_range[1]
	for i,(c,v) in enumerate(zip(coords,values)):
		if np.isnan(v[1]):
			color=color_missing
		else:
			r = (v[1]-minv)/(maxv-minv)
			color = [r*cmax+(1-r)*cmin for cmin,cmax in zip(color_min,color_max)]
		if np.isnan(v[0]):
			width=1
			length = size_missing
		else:
			width=linewidth
			length=v[0]
		ax.add_line(mpl.lines.Line2D([coords[i,0],coords[i,0]], [coords[i,1],coords[i,1]+length], color=color, zorder=2, linewidth=width, solid_capstyle = 'butt'))

	plt.xlim([-1,2**gen])
	plt.ylim([np.nanmax(coords[len(values)/2:len(values),1]+values[len(values)/2:,0])+size_missing/2,-size_missing/2])
	


# yseries: an array the length of objects that holds the series as integers (like a clustering)
# xbins: 'auto' or int (num of bins) or list of values. Only for numeric attributes.
# density = "all" -> sum of everything is 1
#                 "each" -> sum of each bin is 1 (usefull only if y attribute)
#                 "none" or False
# width= float between 0 and 1. Width ratio of an individual bins in its range.
def histplot(objects, xname, yname=None, xbins='auto', xrange='auto', ybins='auto', yrange='auto', yseries=None, stacked=False, density="no", width = 0.8, print_numbers=False, legend=False, newfigure=False):
	
	"""TODO"""
	
	if newfigure:
		plt.figure("{}{}".format(xname," by {}".format(yname) if yname is not None else ""))
	
	# Retrieving values
	_xvalues=data.table(objects, [xname])[:,0]
	if yseries is not None:
		_yvalues = yseries
		yname="subgroup"
	elif yname is not None:
		_yvalues=data.table(objects, [yname])[:,0]
	
	# Removing nones
	nonemask=np.array([v is not None for v in _xvalues],dtype=bool)
	if yname is not None:
		nonemask = nonemask*np.array([s is not None for s in _yvalues], dtype=bool)
		_yvalues=_yvalues[nonemask]
	_xvalues=_xvalues[nonemask]
	
	#Case when attribute y has too many values
	if yname is not None and ybins=='auto':
		_u_yvalues=np.unique(_yvalues)
		if len(_u_yvalues)>=10 and data._get_attribute_type(yname,objects[0])=="nominal":
			raise ExceptException("Too many distinct values for attribute {}".format(yname))
	
	#Building y-series.
	if yname is not None:
	
		_series = {}
		
		if yname=="subgroup" or data._get_attribute_type(yname,objects[0])=="nominal":
			for u in _u_yvalues:
				_series[u] = _xvalues[_yvalues==u]
			
		elif data._get_attribute_type(yname,objects[0])=="numeric":
			if yrange == 'auto':
				yrange = [np.min(_yvalues),np.max(_yvalues)]
			if ybins == 'auto':
				ybins = 10
			if isinstance( ybins, (int,long)):
				ybins = np.arange(0,ybins+1,dtype=float)*float(yrange[1]-yrange[0])/(ybins) + yrange[0]
			for low,high in zip(ybins[:-2],ybins[1:-1]):
				_series["[{:.4f} ; {:.4f}[".format(low,high)] = _xvalues[(_yvalues>=low)*(_yvalues<high)]
			_series["[{:.4f} ; {:.4f}]".format(ybins[-2],ybins[-1])] = _xvalues[(_yvalues>=low)*(_yvalues<=high)]
			
			
	#Counting histograms for numeric x
	if data._get_attribute_type(xname,objects[0])=="numeric":
	
		#Auto-compute range for numeric attributes
		if  xrange=='auto':
			xrange = [np.min(_xvalues),np.max(_xvalues)]
		#When xbins are not provided, does uniform binning:
		#Generates $xbins$ bins of same size in the given range
		if xbins == 'auto':
			xbins = 10
		if isinstance( xbins, (int,long)):
			xbins=np.arange(0,xbins+1,dtype=float)*float(xrange[1]-xrange[0])/(xbins) + xrange[0]
			
		#Building histograms
		if yname is None:
			_hist=np.array(np.histogram(_xvalues, bins=xbins)[0], int)
		else:
			_hist=np.array([np.histogram(_series[v], bins=xbins)[0] for v in _series], int)
			
	
	#Counting histograms for nominal x
	elif data._get_attribute_type(xname,objects[0])=="nominal":
		
		_u_xvalues = np.unique(_xvalues)
		#cancel if too many values?
		
		#Building histograms
		if yname is None:
			_hist=np.array([np.sum(_xvalues==x) for x in _u_xvalues], int)
		else:
			_hist=np.array([[np.sum(_series[v]==x) for x in _u_xvalues] for v in _series], int)
			
		
		xbins=np.arange(len(_u_xvalues)+1)
		plt.xticks(xbins[:-1]+0.5, _u_xvalues)
		
		
	#Normalisation
	if density == "all" and yname is not None:
		_bars = _hist.astype(float)/np.reshape(np.sum(_hist, axis=1),(-1,1))/_hist.shape[0]
	elif (density == "each" and yname is not None) or (density == "all" and yname is None):
		_bars = _hist.astype(float)/np.reshape(np.sum(_hist, axis=0),(1,-1))
	else:
		_bars=_hist
	
	#Colors
	colors="brgycmk"
	if yname is not None:
		while(len(colors)<len(ybins)-1):
			colors=colors+colors
		colors=colors[:len(ybins)-1]
	
	#Plotting
	plt.xlabel(xname)
	for b in range(len(xbins)-1):
		_bottom = 0
		_width = float(width*(xbins[b+1]-xbins[b]))
		if (not stacked) and yname is not None:
			_width = _width/(len(ybins)-1)
		_left = xbins[b] + (xbins[b+1] - xbins[b])*(1 - width)/2
		
		#No series
		if yname is None:
			# _height = _bars[0,b]
			_height = _bars[b]
			if not isnan(_height) and _height>0:  #NaN possible?
				plt.bar(_left,_height,_width,_bottom, align='edge', color=colors[0])
				if print_numbers:
					plt.text( _left+_width/2., _bottom + _height/5, str(_hist[b]), color='w', ha="center")
		else:
			for i,s in enumerate(_series):
				_height = _bars[i,b]
				plt.bar(_left,_height,_width,_bottom, align='edge', color=colors[i])
				if not isnan(_height) and _height>0:
					if print_numbers:
						plt.text( _left+_width/2., _bottom + _height/5, str(_hist[i,b]), color='w', ha="center" )
				if stacked and not isnan(_height):
					_bottom=_bottom+_height
				else:
					_left = _left+_width
						
		
	if yname is not None:
		plt.legend(_series, title=yname)
		
	return plt.gca()
	
		
def scatterplot(cells, xname, yname, mask=None, color_attr=None,color_mask=None, regression=None):
	
	"""TODO"""

	if mask is None:
		mask=np.ones(len(cells),dtype='bool')
	if color_mask is None:
		color_mask=np.zeros(len(mask),dtype='int')
	if color_attr is not None:
		#numerical
		colors = data.table(cells[mask], [color_attr])
		color_bins=np.histogram(colors, bins = min(10, len(np.unique(colors))-1))[1]
		color_mask = (colors>color_bins).sum(axis=1)
		#categorical
		#TODO
		
	_data=data.table(cells[mask], [xname,yname])
	submask=np.array([d[0] is not None and d[1] is not None for d in _data])
	X=np.array(_data[submask,0], dtype='float')
	Y=np.array(_data[submask,1], dtype='float')
	color_mask=color_mask[submask]
	
	labels=np.unique(color_mask)
	
	for l in labels:
		plt.plot(X[color_mask==l],Y[color_mask==l], '.')
	plt.xlabel(xname)
	plt.ylabel(yname)
	
	if regression=="linear":
		clf=linear_model.LinearRegression()
		clf.fit([[x] for x in X],Y)
		minmax=[[X.min()],[X.max()]]
		plt.plot(minmax, clf.predict(minmax), 'k-')
		plt.figtext(0.2,0.95,"%s = %0.3g + %0.3g x %s  -  R2=%0.3g"%(yname, clf.intercept_, clf.coef_[0], xname, clf.score([[x] for x in X],Y)))

def boxplot(objects, xname, yname):
	
	"""TODO"""
	
	_data=data.table(objects,[xname, yname])
	_mask=np.array([d[0] is not None and d[1] is not None for d in _data])
	_data=_data[_mask]
	_boxes=[_data[np.array([x[0]==a for x in _data]),1] for a in np.sort(np.unique(_data[:,0]))]

	# plt.figure()
	plt.boxplot(_boxes, positions=np.sort(np.unique(_data[:,0])).astype('int'))
	plt.xlabel(xname)
	plt.ylabel(yname)

