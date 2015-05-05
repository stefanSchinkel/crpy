import crpy
X = sin(linspace(1,6*pi,500))
X = loadtxt('data.dat')
rp = crpy.crp(X,[],3,2,.1)
crpy.crqa(X,[],3,2,.1)

X = sin(linspace(1,6*pi,10))
rp = crpy.crp(X,[],3,2,.1)
crpy.crqa(X,[],3,2,.1)

