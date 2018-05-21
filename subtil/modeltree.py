import numpy as np
import time
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.tree._tree import Tree
from sklearn.tree import DecisionTreeRegressor

#for test
X=np.array([[1.,1.,1.4,1.9,2.1,2.9,2.9,3.8,3.9,4.5,5.,5.6],[1.1,2.,2.2,1.,2.3,0.9,2.5,1.2,1.8,2.,1.1,1.]]).T
y=np.array([1.,2.,2.3,1.7,3.4,2.3,4.5,3.4,5.4,6.7,3.6,4.5])

class ModelTreeRegressor(DecisionTreeRegressor):

	"""TODO"""
	
	def __init__(self, max_depth=3, min_samples_leaf=25, debug=False):
		self.max_depth=max_depth
		self.min_samples_leaf=min_samples_leaf
		self.tree=Tree(np.size(X,axis=1),np.array([1]),1)
		self.debug=debug
		
	def fit(self,X,y):
		
		self.children_left=[-1]
		self.children_right=[-1]
		self.node_count=1
		self.feature=[]
		self.threshold=[]
		self.n_node_samples=[]
		self.linear_models=[Ridge(alpha=0.01)]
		self.linear_models[0].fit(X,y)
		self.mse = [((self.linear_models[0].predict(X)-y)**2).sum()]
		
		self._split(X,y,0)
		state={'node_count':self.node_count,
				'values':np.array([[self.mse]],order='C').reshape((-1,1,1)),
				'nodes':np.array([(a,z,e,r,t,y,u) for a,z,e,r,t,y,u in zip(self.children_left,self.children_right,self.feature,self.threshold,self.mse,self.n_node_samples,self.n_node_samples)],dtype=self.tree.__getstate__()['nodes'].dtype)}
		# return state
		self.tree.__setstate__(state)
		return self
		
			
	def _split(self, X, y, node, depth=0):
		
		if depth>=self.max_depth:
			self.feature.append(-1)
			self.threshold.append(0.)
			self.n_node_samples.append(len(X))
			return
		
		if self.debug:
			print("")
			print("left {}".format(self.children_left))
			print("right {}".format(self.children_right))
			print("splitting node {}, {} points, mse={}".format(node, len(X),self.mse[node]))
			# print X
			# print y
		
		best_a=-1
		best_threshold=0
		best_mse=self.mse[node]
		best_mask=None
		for a in range(np.size(X,axis=1)):
			if self.debug: print("attribute {}".format(a))
			arg=np.argsort(X[:,a])
			mask=np.ones(len(X),dtype=bool) #1 to the right, 0 to the left
			for i in range(0,len(X)-self.min_samples_leaf):
				if i<self.min_samples_leaf or X[arg[i],a]==X[arg[i+1],a]:
					mask[arg[i]]=False
					continue
				# if self.debug: print("  threshold {}".format(X[arg[i],a]))
				# L1=LinearRegression()
				# L2=LinearRegression()
				L1=Ridge(alpha=0.01)
				L2=Ridge(alpha=0.01)
				L1.fit(X[mask==False],y[mask==False])
				L2.fit(X[mask],y[mask])
				y1=L1.predict(X[mask==False])
				y2=L2.predict(X[mask])
				mse1=((y1-y[mask==False])**2).sum()
				mse2=((y2-y[mask])**2).sum()
				
				# print "    {} points, mse1={}".format(i,mse1)
				# print "    {} points, mse2={}".format(len(X)-i,mse2)
				mse=(float(i)/len(X))*(((y1-y[mask==False])**2).sum())+(float(len(X)-i)/len(X))*(((y2-y[mask])**2).sum())
				# if self.debug: print("  total mse={}".format(mse))
				if mse<best_mse:
					# print mask
					best_a=a
					best_threshold=(X[arg[i],a]+X[arg[i+1],a])/2
					best_mse=mse
					best_mse1=mse1
					best_mse2=mse2
					best_l1=L1
					best_l2=L2
					best_mask=np.array(mask)
				mask[arg[i]]=False
				
				####### ?????????????? To remove?########
				if self.debug: time.sleep(0.001)
				
		if best_a==-1:
			self.feature.append(-2)
			self.threshold.append(-2.)
			self.n_node_samples.append(len(X))
			return
		
		self.feature.append(best_a)
		self.threshold.append(best_threshold)
		self.n_node_samples.append(len(X))
		
		#create left children node
		self.children_left.append(-1)
		self.children_right.append(-1)
		self.mse.append(best_mse1)
		self.linear_models.append(best_l1)
		self.children_left[node]=self.node_count
		self.node_count=self.node_count+1
		self._split(X[best_mask==False],y[best_mask==False],self.node_count-1,depth=depth+1)
		
		#create left children node
		self.children_left.append(-1)
		self.children_right.append(-1)
		self.mse.append(best_mse2)
		self.linear_models.append(best_l2)
		self.children_right[node]=self.node_count
		self.node_count=self.node_count+1
		self._split(X[best_mask],y[best_mask],self.node_count-1,depth=depth+1)
		
		
	def predict(self,X):
		if len(np.shape(X))==1:
			X=np.array(X).reshape(1,-1)
			
		predicted=np.empty(len(X))
		nodes=self.tree.apply(np.array(X,dtype=np.float32))
		# print nodes
		for i,n in enumerate(nodes):
			predict=self.linear_models[n].predict(X[i])
			if predict<0:
				predict=0
			elif predict>1:
				predict=1
			predicted[i]=predict
			
		return predicted
		

