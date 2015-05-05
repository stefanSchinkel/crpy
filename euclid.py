def euclidDistance(point1,point2):

	dist = 0
	for n in range(len(element1)):
		dist = dist + pow((element1[n] - element2[n]),2)
	#print 'dist' + dist

dim = 3;tau=3;
X = sin(linspace(1,6*pi,20))
Xemb = embedData(X,dim,tau)

for i in range(size(Xemb,0)): print Xemb[i,:];
	distance = 0 	
	for j in range(dim):
	
