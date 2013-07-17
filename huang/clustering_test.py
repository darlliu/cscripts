from clustering import svdpca,kmeans
import os, sys, string, numpy

ss=[[1,2],[4,2],[300,200]]
ll=svdpca(str(ss), 3, 2);
print ll;
ss=[[300,200],[1,2],[4,2]]
ll=kmeans(str(ss), 2, 2);
print ll;
#exchange
#sss = [[1,2],[4,2],[3,1]]
#ll=svdpca(str(ss), 3,2);
#print ll;
#ll=svdpca(str(sss), 3,2);
#print ll;
