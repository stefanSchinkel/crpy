#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
	Purpose:
	
	This class provides a set of functions necessary for Recurrence 
	Quantification Analysis. i.e:
		- computation of distance matrices & recurrence plots
		- quantification of recurrence plot structure (not yet, sorry)
	
	Requirements:
		- numpy  % for MATLAB-like operations (http://numpy.scipy.org/)
		- pylab  % for plotting stuff (http://matplotlib.sourceforge.net/)
	
	Notes:
		Core functions (lowercase):
			- should be called directly
		Helper functions (camelCase):
			- get called by Core funtion set. 
			-> should NOT be called on their own
		Plotting functions (show*):
			- will probably get a special wrapper function

	Implemented Features: 		
		Core:		
			- crp(X, Y, dim, tau, eps, norm='max', normalise=True)
		Helper: 
			- normaliseData(X)
			- embedData(X, dim, tau)
 			- makeDistMatrix(X, Y, dim, tau)
 		Plotting: 
 			- showDistMatrix(distMatrix)
 			- showRP(RP)
 	ToDo:
 		Core:
	 		- crqa(X, Y, dim, tau, eps, norm='max', normalise=True)
			- rqaci(X, Y, dim, tau, eps, norm='max', normalise=True)
			- rqabounds(X, Y, dim, tau, eps, norm='max', normalise=True)
		Helper: 
			- applyTheiler(RP)
	 		- getVertLines(RP)
 			- getDiagLines(RP)
 			- quantifyLineDistribution(pLines) -> should be same for diag/vert
	 		- bootstrapRQAMeasures()
	 	Plotting:
	 		- plotRQAMeasures()
 		
	MENTAL NOTES: 
		Maybe we should loose makeDistMatrix. 
		It is kind of redundant, but may help 
		to save A LOT of memory.
			
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
##########################
##  	Imports 		##
##########################

try:
	from numpy import *
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

#the external functions
from helper import *
from plots import * 
from other import *


##################################
##  	Core Functions			##
##################################

def crp(X,Y,dim,tau,eps,norm='max',normalise=True):		
	"""
	Computes a recurrence plot with the given embedding parameters,
	threshold and norm and returns it to the caller. For the sake 
	of memory saving no distance matrix is created. The data is embedded
	within the function. Multi-column input not yet supported. 
	"""	
	
	# RP/CRP switch	
	if len(Y) == 0: Y = X

	# for alloc and looping
	lenData = len(X)	
	lenEmbData = lenData - (dim-1)*tau		

	#normalise
	if normalise: 
		X = normaliseData(X)
		Y = normaliseData(Y)

	#embed
	Xemb = embedData(X,dim,tau)
	Yemb = embedData(Y,dim,tau)

	#alloc memory for CRP matrix
	rp = zeros( (lenEmbData,lenEmbData), dtype=int8)	
	
	#loop
	for i in range(lenEmbData):
		
		
		xDiff = abs(Xemb[i,:] - Yemb)
		# "switch norms"
		if norm == 'max':	

			
			rp[xDiff.max(axis=1) < eps,i] = 1
			#diff = abs(Yemb) - matlib.repmat(abs(Xemb[i]),lenEmbData,1)
			#rp[:,i] = diff.max(axis = 1) < eps

		elif norm == 'min':

			diff = Yemb - matlib.repmat(Xemb[i],lenEmb,1)				
			rp[:,i] = abs(diff.min(axis = 1)) < eps
		
		else :				
			raise Exception,"Only maximum and minimum norm supported now"
			
	#return to caller
	return rp

def crqa(X,Y,dim,tau,eps,norm='max',theiler=1,lMin=2,vMin=2,normalise=True):
	"""
	This function computes an (C)RP from the given data
	and quantifies it. The complexity measures a returned
	to the caller in rqa:
	
		rqas[0] == RR
		rqas[1] == DET*
		rqas[2]	== L
		rqas[3] == Lmax
		rqas[4] == ENT*
		rqas[5] == LAM*	
		rqas[6] == TT
		rqas[7] == Vmax

		* == not implemented yet
	
	windowing is not supported yet. 
	"""
	
	# RP/CRP switch	
	if len(Y) == 0: Y = X
	
	rp = crp(X,Y,dim,tau,eps,norm,normalise)
	
	rqa = qualifyRP(rp,theiler,lMin,vMin)

	print "RR\tDET\tL\tLmax\tENT\tLam\tTT\tVmax\t"
	print "%3.4f\t%3.4f\t%3.4f\t%3.1f\t%3.4f\t%3.4f\t%3.4f\t%3.1f\t" % \
	 	(rqa[0],rqa[1],rqa[2],rqa[3],rqa[4],rqa[5],rqa[6],rqa[7])


			
##################################
##  		Sample Code			##
##################################

if __name__ == "__main__":


	from helper import *
	from plots import * 
	from other import *
	print dir()
	
	print "\tHey there, looks like I'm run standalone. "
	print "\tI'll show you a little demo then.\n";
	
	X = sin(linspace(1,6*pi,500))
	Y = sin(linspace(1,6*pi,500))
	dMat = makeDistMatrix(X,[],3,2)
	showDistMatrix(dMat)	
	rp = crp(X,Y,3,2,.1)
	showRP(rp)
	
	print "demo is done"
