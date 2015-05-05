#!/usr/bin/env python
# -*- coding: utf-8 -*-
from numpy import *
from numpy import matlib
import pylab

"""
	Purpose:
	
	These are the plotting functions for the crpy module
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


##################################
##  	Plotting Functions		##
##################################

def showDistMatrix(distMatrix):		
	pylab.matshow(distMatrix)
 	ax = pylab.gca() 
 	bottom, top = ax.get_ylim()
	ax.set_ylim(top, bottom)		
	pylab.colorbar() 
	pylab.show()


def showRP(RP):

	pylab.matshow(- RP)
 	ax = pylab.gca() 
 	bottom, top = ax.get_ylim()
	ax.set_ylim(top, bottom)
	pylab.axis('tight')
	pylab.gray()
	pylab.show()
