#how to get the VERTICAL RP lines:
from numpy import *
# pseudo RP 
RP = eye( 5 )

#pad with zeros before vectorising
padding = zeros( (5,1), dtype=int8)
RP = concatenate( (padding,RP,padding), axis=1)

#concat for maxConsElements
lines  = RP.reshape( (RP.size) )  

#maxConsElements
tmp = diff(lines)
ind1 = array( where(tmp == 1), dtype = int32 )
ind2 = array( where(tmp == -1), dtype = int32 )
vertLines = ind2-ind1

print vertLines

vec = array ([0,1,0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0 ])
tmp = diff(vec,axis=0)
ind1 = array( where(tmp == 1), dtype = int )
ind2 = array( where(tmp == -1), dtype = int )
distLines = ind2-ind1
