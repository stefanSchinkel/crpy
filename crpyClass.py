#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
	Purpose:
	
	This class provides a set of functions necessary for Recurrence 
	Quantification Analysis. I.e:
		- normalising data
		- embedding of timeseries
		- computation of distance matrices & recurrence plots
		- quantification of recurrence plot structure (not yet, sorry)
	
	Requirements:
		- numpy  % for MATLAB-like operations (http://numpy.scipy.org/)
		- pylab  % for plotting stuff (http://matplotlib.sourceforge.net/)
	
	Notes:
		lowercase functions are convenience function : 
			they are the only ones meant to be called directly	
		camelCase functions are helper functions:
			they are not meant to run indepently, the can be 
			called though - but be sure you know what you're doing
			
	* Copyright (c) 2008, Stefan Schinkel, University of Potsdam
	* All rights reserved.
	*
	* Redistribution and use in source and binary forms, with or without
	* modification, are permitted provided that the following conditions are met:
	*     * Redistributions of source code must retain the above copyright
	*       notice, this list of conditions and the following disclaimer.
	*     * Redistributions in binary form must reproduce the above copyright
	*       notice, this list of conditions and the following disclaimer in the
	*       documentation and/or other materials provided with the distribution.
	*     * Neither the name of the University of Potsdam nor the
	*       names of its contributors may be used to endorse or promote products
	*       derived from this software without specific prior written permission.
	*
	* THIS SOFTWARE IS PROVIDED BY Stefan Schinkel ''AS IS'' AND ANY
	* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
	* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	* DISCLAIMED. IN NO EVENT SHALL Stefan Schinkel BE LIABLE FOR ANY
	* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
	* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
	* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

	$Log$


"""

class crpy:
	""" 
	This is the crpy class. No, its not really a module. 
	Therefore please import it as: from crpy import *.
	The class provides a basic set of functions for computing
	and displaying RPs and CRPs.
	"""
	
		
	def __init__(self,X,Y=[]):	
		"""
		Initialises all params and flags.
		Stores X and Y in object.
		"""
		
		__name__ = 'crpy'
		# init non-boolean helpers
		self.dim = []
		self.tau = []
		self.eps = []
		self.norm = ''
		
		# all boolean helpers, for later use
		self.flagIsCRP = bool()
		self.flagIsNormalised = False
		self.flagIsEmbedded = False
		self.flagIsDistMatrix = False
		self.flagIsRP = False
		
		# assign data
		self.X = X;
		
		if len(Y) == 0:
			self.Y = X
			# normalise/embed only once
			self.flapIsCRP = False;

		else: 
			self.Y = Y
			# normalise/embed only twice
			self.flapIsCRP = True;			
			if len(self.X) != len(self.Y):
				raise Exception,"X and Y must be of the same length"	

	
	def normaliseData(self):
		"""
		Normalises X and Y to mu=0 & sigma = 1	(in-place)
		"""
		
		# mu & sigma
		meanX = self.X.mean()
		stdX = self.X.std()

		if not self.flagIsCRP:
		
			self.X -= meanX;self.X /= stdX;
			self.Y = self.X

		else:
		
			meanY = self.Y.mean()
			stdY = self.Y.std()
			
			self.X -= meanX;self.X /= stdX
			self.Y -= meanY;self.Y /= stdY
		

	def embedData(self,dim,tau):
		"""
		Embeds the vectors a X and Y with given 
		embedding dimension and delay
		"""
		
		# data length		
		self.lenData = len(self.X)

		# consistency check 		
		if (dim-1)*tau >= self.lenData:
			raise ValueError,"Data not long enough for embedding parameters"		

		self.lenEmbData = self.lenData - (dim-1)*tau		

		# store for later plotting/debugging
		self.dim = dim;
		self.tau = tau;

		if not self.flapIsCRP:
		# only embedded once since X == Y
			
			# alloc memory		
			self.Xemb = zeros( (self.lenData - (dim-1)*tau, dim), dtype=float32 )

			# assign values
			for i in range(dim):

				self.Xemb[:,i] = self.X[i*tau : self.lenData -( (dim-1)-i ) *tau]

			self.Yemb = self.Xemb

		else:
			# alloc memory		
			self.Xemb = zeros( (self.lenData - (dim-1)*tau, dim), dtype=float32 )
			self.Yemb = zeros( (self.lenData - (dim-1)*tau, dim), dtype=float32 )

			# assign values
			for i in range(dim):

				self.Xemb[:,i] = self.X[i*tau : self.lenData -( (dim-1)-i ) *tau]
				self.Yemb[:,i] = self.Y[i*tau : self.lenData -( (dim-1)-i ) *tau]
		
	def makeDistMatrix(self,dim,tau,norm='max',normalise=True):
		"""
		Computes the chosen distance between all
		points in X and Y and stores the matrix
		in the class.
		"""	

		#normalise
		if normalise: self.normaliseData()
		#embed
		self.embedData(dim,tau)

		# alloc output matrix
		self.distMatrix = empty( (self.lenEmbData,self.lenEmbData), dtype=float32)	

		for i in range(self.lenEmbData):

			# "switch norms"
			if norm == 'max':	
			
				diff = self.Yemb - matlib.repmat(self.Xemb[i],self.lenEmbData,1)				
				self.distMatrix[:,i] = 	abs(diff.max(axis = 1))		
				
			elif norm == 'min':

				diff = self.Yemb - matlib.repmat(self.Xemb[i],self.lenEmbData,1)			
				self.distMatrix[:,i] = 	abs(diff.min(axis = 1))

			else:

				raise Exception,"Only maximum and minimum norm supported now"
	
		# for short-circuiting				
		self.flagIsDistanceMatrix = True	


	def crp(self,dim,tau,eps,norm='max',normalise=True):		
		"""
		Computes a recurrence plot with the given embedding parameters,
		threshold and norm. 
		"""	
		#normalise
		if normalise: self.normaliseData()
		#embed
		self.embedData(dim,tau)

		#and now for the crp
		#alloc memory for CRP matrix
		self.rp = 	zeros( (self.lenEmbData,self.lenEmbData), dtype=int8)	
		
		#loop
		for i in range(self.lenEmbData):

			# "switch norms"
			if norm == 'max':	
			
				diff = self.Yemb - matlib.repmat(self.Xemb[i],self.lenEmbData,1)
				self.rp[:,i] = 	abs(diff.max(axis = 1)) < eps

			elif norm == 'min':

				diff = self.Yemb - matlib.repmat(self.Xemb[i],self.lenEmb,1)				
				self.rp[:,i] = 	abs(diff.min(axis = 1)) < eps
			
			else :				
				raise Exception,"Only maximum and minimum norm supported now"
				
		# for short-circuiting		
		self.flagIsRP = True			

	def makeOPRP(self):		
		
		# DO HERE
		# symbolise timeseries
		# match symbols for OPRP
		pass


##################################
##  Data handling/ Class I/O	##
##################################

	def getDistMatrix(self):
		"""
		Return distance matrix to caller.
		"""	
	
		if not self.flagIsDistMarix:
	
			print "No Distance Matrix computed yet. Sorry."
			return

		else:
			return self.DistMatrix
		
	def getRP(self):
		"""
		Return  matrix to caller.
		"""	
	
		if not self.flagIsRP:
			print "No RP computed yet. Sorry."
			return
		else:
			return self.RP
	

##################################
##  	Plotting Functions		##
##################################

	def showDistMatrix(self):		
		pylab.matshow(self.distMatrix)
	 	ax = pylab.gca() 
	 	bottom, top = ax.get_ylim()
		ax.set_ylim(top, bottom)		
		pylab.colorbar() 
		pylab.show()


	def showRP(self):

		pylab.matshow(- self.rp)
	 	ax = pylab.gca() 
	 	bottom, top = ax.get_ylim()
		ax.set_ylim(top, bottom)
		pylab.axis('tight')
		pylab.gray()
		pylab.show()


##################################
##  		Sample Code			##
##################################


if __name__ == "__main__":

	# default import
	import sys,time

	print "\tHey there, looks like I'm run standalone. "
	print "\tI'll show you a little demo then.\n";

	try:
		from  numpy import *
	except ImportError:
		error('Cannot load numpy module!')
	try:
		from  numpy import matlib
	except ImportError:
		error('Cannot load numpy-module matlib!')
	try:
		import pylab
	except ImportError:
		error('Cannot load pylab module!')
	
	# for illustration a sine wave	
	print "\tI'll prep some sample data. A sine X & a cosine Y.\n" 
	X = sin(linspace(1,6*pi,500))
	Y = sin(linspace(1,6*pi,500))
	
	print "\tpython>crp = crpy(X,Y) instatiates an crpy Object\n"
	crp = crpy(X,Y)

	print "\tpython>crp.makeDistMatrix(3,2)"
	print "\twill compute the distance matrix with dim=3,tau=2 (maximum norm)."
	print "\tThe data is normalised by default.\n"
	crp.makeDistMatrix(3,2)

	print "\tLets look at the Matrix....\n\tClose the Figure to proceed!\n"	
	crp.showDistMatrix()	

	print "\tpython>crp.crp(3,2,.1)"
	print "\twill compute a maximum norm CRP with dim=3,tau=2 and eps=.1"
	print "\tThe data is normalised by default.\n"

	crp.crp(3,2,.1);	

	print "\tLets look at the RP....\n\tClose the Figure to proceed!\n"
	crp.showRP()

	print "\tThat's all so far. More to come."
	
	