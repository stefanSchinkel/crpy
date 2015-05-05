#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
	Purpose:
	
	These are the helper functions for the crpy module
	See crpy for details. 

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
from numpy import *
from numpy import matlib
import pylab

##################################
##  	Helper Functions		##
##################################


def normaliseData(X):
	"""
	Normalises X and Y to mu=0 & sigma = 1	(in-place)	
	"""
	
	# mu & sigma
	meanX = X.mean()
	stdX = X.std()

	X -= meanX;
	X /= stdX;

	return X

def embedData(X,dim,tau):
	"""
	Embeds vector a X with given 
	embedding dimension and delay
	and returns it to caller
	"""
	
	# data length		
	lenData = len(X)

	# consistency check 		
	if (dim-1)*tau >= lenData:
		raise ValueError,"Data not long enough for embedding parameters"		

	lenEmbData = lenData - (dim-1)*tau		
		
	# alloc memory
	Xemb = zeros( (lenData - (dim-1)*tau, dim), dtype=float32 )

	# assign values
	for i in range(dim):

		Xemb[:,i] = X[i*tau : lenData -( (dim-1)-i ) *tau]

	#return to caller
	return Xemb

def getDiagLines(RP):
	"""
	Extract the histogramme of all diagonal lines
	in a recurrence plot. Note: numpy does not 
	have convenient spdiags(), therefore the
	crappy looping.
	"""
	
	# no. of diagonals on one axis
	nDiags = alen(RP);
	
	#alloc empty array, force int
	lines = array([], dtype=int32)

	# get all diagonals and extract line lengths	
	for offset in range(-nDiags +1,nDiags):
		lines = append( lines, maxConsElements(diag(RP,offset)) )
		
	return lines

def getVertLines(RP):
	"""
	Extract the histogramme of all vertical lines
	in a recurrence plot.
	"""

	#alloc empty array, force int
	lines = array([], dtype=int32)
	
	#extract all columns individually
	for i in range(alen(RP)):
		lines = append( lines, maxConsElements(RP[:,i]) ) 

	return lines

def maxConsElements(vector):
	"""
	Returns the length of all consecutive elements in vector. 
	"""
	
	#vector MUST start/end with 0	
	vector = concatenate( ([0],vector,[0]) )
	
	tmp = diff(vector)
	ind1 = array( where(tmp == 1), dtype = int32 )
	ind2 = array( where(tmp == -1), dtype = int32 )
	lines = ind2-ind1	
	
	return lines

def qualifyRP(RP,theiler=1,lMin=2,vMin=2):
	"""
	This functions takes an RP matrix as an input
	and computes the complexity measures
	"""
	
	#array to hold results
	rqa = zeros(8)
		
	#theiler exclusion here please
	if theiler > 0:
		RP = triu(RP,theiler)+tril(RP,-theiler)
	
	# all recurrent points	
	allPoints = sum(RP, dtype=float32)

	# get the line distibution	
	diagLines = getDiagLines(RP)
	vertLines = getVertLines(RP)

	# and exclude lines that are too short
	diagLines = diagLines[diagLines >= lMin]
	vertLines = vertLines[vertLines >= vMin]
	
	# compute the complexity measures	
	
	# RR
	rqa[0] = allPoints/RP.size
	

	#DET 
	rqa[1] = sum(diagLines) / allPoints	
	
	# L und Lmax
	rqa[2] = diagLines.mean()
	rqa[3] = diagLines.max()
	
	#ENT
	(p,h) =	histogram(diagLines,new=True)
	rqa[4] = shannon(p/allPoints)
	
	#LAM 
 	rqa[5] = sum(vertLines) / allPoints

	# TT and Vmax
	rqa[6] = vertLines.mean()
	rqa[7] = vertLines.max()

	
	return rqa
	
def shannon(p):
 	"""
	Compute Shannon entropy of random variable with probability p
	"""
	#log 0 == -inf
	p = p[p > 0]
	
	return sum(-p*log2(p))




