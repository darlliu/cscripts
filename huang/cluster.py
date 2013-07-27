from clustering import svdpca, kmeans
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os, sys, string, math
import numpy as np
import cPickle as pkl
import crick.viz.js.google_chart_js as chart
import matplotlib.cm as cm
from wrapper import convertIDs
import crick.viz.js.google_chart_js as chart
import crick.annotate.cybert.crick_cybert_annotator as anno
import crick.builders.protein_protein_edges.StringBuilder as sb
from crick.utils.tabular import tabular
from crick.const import species

targets = os.listdir("./csv")
#print targets
tags = open("./csv/mean_all.txt","r").read().split("\n")
tags = sorted([i.strip() for i in tags])
ids = convertIDs (fromType="UNIPROT_ACCESSION", toType="OFFICIAL_GENE_SYMBOL",ids=list(tags))
IDS = dict(ids)
print len(tags)
htags= [IDS[i] if IDS[i]!="N/A" else i for i in tags]
open("./sources.list","w").write("\n".join(htags+["Identified Targets"]))
print len(tags)
#print tags

PVALS=pkl.load(open("./pkl/NatDatapvals_all_all.pkl","rb"))
MEANS=pkl.load(open("./pkl/NatDatameans_all_all.pkl","rb"))
CMEANS=pkl.load(open('./pkl/NatDatacmeans_all_all.pkl',"rb"))
#print CMEANS
NMEANS=pkl.load(open('./pkl/NatDatanmeans_all_all.pkl',"rb"))
BAITS=["Pre1","Pre10","Rpt1","Rpn1","Rpt5","Rpt6","Rpn11"]+["Base","Lid","20S"]

nat_pfxs= ["Nat_"+i for i in BAITS]
den_pfxs= ["Den_"+i for i in BAITS]
#nat_gpfxs= ["Nat_"+i for i in GBAITS]
#den_gpfxs= ["Den_"+i for i in GBAITS]

#from here everything is native
all_means = [[float(MEANS[(tag,pfx)]) if (tag,pfx) in MEANS else 0.0 for pfx in nat_pfxs] for tag in tags]
all_means = tabular (nat_pfxs, all_means)
all_cmeans = [[float(CMEANS[(tag,pfx)]) if (tag,pfx) in CMEANS else 0.0 for pfx in BAITS] for tag in tags]
all_cmeans = tabular (nat_pfxs,all_cmeans )
all_nmeans = [[float(NMEANS[(tag,pfx)]) if (tag,pfx) in NMEANS else 0.0 for pfx in BAITS] for tag in tags]
all_nmeans = tabular(nat_pfxs, all_nmeans)
all_pvals = [[float(PVALS[(tag,pfx)]) if (tag,pfx) in PVALS else 1.0 for pfx in nat_pfxs] for tag in tags]
all_pvals = tabular (nat_pfxs, all_pvals)
#print all_pvals
#print all_pvals
sig_flags = [[float(i)<0.05 for i in j] for j in all_pvals]
#print sig_flags

filtered_means=[]
filtered_pvals=[]
#now produce filtered data
SUBNETS={}


for mean_row, pval_row, flag_row in zip (all_means, all_pvals, sig_flags):
    mm= []
    pp= []
    for m, p , f in zip(mean_row, pval_row, flag_row):
        if f:
            mm.append(float(m))
            pp.append(float(p))
        else:
            mm.append(0.0)
            pp.append(1.0)

    filtered_means.append(mm)
    filtered_pvals.append(pp)

print filtered_means==all_means
print filtered_pvals==all_pvals
#now produce the reduced data

# now split columns

baits_means = [i[:-3] for i in all_means]
baits_cmeans = [i[:-3] for i in all_cmeans]
baits_nmeans = [i[:-3] for i in all_nmeans]
baits_pvals = [i[:-3] for i in all_pvals]
baits_filtered_means = [i[:-3] for i in filtered_means]
baits_filtered_pvals = [i[:-3] for i in filtered_pvals]
print nat_pfxs[:-3], nat_pfxs[-3:]

grouped_means = [i[-3:] for i in all_means]
grouped_pvals = [i[-3:] for i in all_pvals]
grouped_filtered_means = [i[-3:] for i in filtered_means]
grouped_filtered_pvals = [i[-3:] for i in filtered_pvals]

def do_heat (tags, pfxs, data, title="Heatmap"):
    # assemble the rows based on pfxs
    ch=chart.proper_heat(title)
    ch.load(pfxs, data, tags)
    return

def f(col):
    return [[-math.log(i,10) for i in j] for j in col]

def routine1():
    #plot various heatmaps
    do_heat(htags, nat_pfxs[:-3], f(baits_pvals), "p Values of the baits, raw values, negative log scale")
    do_heat(htags, nat_pfxs[:-3], f(baits_filtered_pvals), "p Values of the baits, filtered values, negative log scale")

    do_heat(htags, nat_pfxs[-3:], f(grouped_pvals), "p Values of the groups, raw values, negative log scale")
    do_heat(htags, nat_pfxs[-3:], f(grouped_filtered_pvals), "p Values of the groups, filtered values, negative log scale")
    do_heat(htags, nat_pfxs[:-3], baits_means, "Means of the baits, raw values")
    do_heat(htags, nat_pfxs[:-3], baits_filtered_means, "Means of the baits, filtered values")
    do_heat(htags, nat_pfxs[:-3], baits_pvals, "p Values of the baits, raw values")
    do_heat(htags, nat_pfxs[:-3], baits_filtered_pvals, "p Values of the baits, filtered values")

    do_heat(htags, nat_pfxs[-3:], grouped_means, "Means of the groups, raw values")
    do_heat(htags, nat_pfxs[-3:], grouped_filtered_means, "Means of the groups, filtered values")
    do_heat(htags, nat_pfxs[-3:], grouped_pvals, "p Values of the groups, raw values")
    do_heat(htags, nat_pfxs[-3:], grouped_filtered_pvals, "p Values of the groups, filtered values")

    do_heat(htags, nat_pfxs, all_means, "Means of the everything, raw values")
    do_heat(htags, nat_pfxs, filtered_means, "Means of everything, filtered values")
    do_heat(htags, nat_pfxs, all_pvals, "p Values of everything, raw values")
    do_heat(htags, nat_pfxs, filtered_pvals, "p Values of everything, filtered values")
    do_heat(htags, nat_pfxs, f(all_pvals), "p Values of everything, raw values, negative log scale")
    do_heat(htags, nat_pfxs, f(filtered_pvals), "p Values of everything, filtered values, negative log scale")
    #print tags[-3], all_means[-2] , filtered_means[-3] ,all_pvals[-3], filtered_pvals[-3]
#routine1()


#kmeans clustering on means and on p values
def k (tags, data, K, tries=5):
    #first get the clustering info
    idxs=kmeans(str(data), K, tries)
    idxs=eval(idxs)
    print len(idxs)==len(tags)
    #for the clusters, get the numbers for eac``
    return idxs

def pca(K, idxs, _data,title, pfxs):
    colors=["r","g","b","k","c","m","y","p"]
    C=colors[:len(pfxs)]
    print C
    #ch=chart.google_line_chart(title="PCA plot of K-Means clusters of {0} with K={1:d}".format(title,K),
            #xlabel="PC1", ylabel="PC2")
    fig=plt.figure(figsize=(12, 5*K/2))
    plt.title("PCA plot of K-Means clusters of {0} with K={1:d}".format(title,K))
    X=[]
    Y=[]
    C=[]
    for i in xrange(K):
        ax=plt.subplot("{0:d}2{1:d}".format(int((K+1)/2),i+1))
        data2=[x[1] for x in filter(lambda x: x[0]==i, zip(idxs, _data))]
        #get data vector
        data=np.array(data2)
        #subtract mean
        mm = np.mean(data)
        mmax = abs(np.amax(data))*1.2
        #mmax=1
        #data = data - mm
        print np.mean(data)
        #get two PCs
        temp_pca=np.array(eval(svdpca(str(data2),len(data2),len(data2[0]))[1]))
        temp_var=np.array(eval(svdpca(str(data2),len(data2),len(data2[0]))[0]))
        print "norms:", np.linalg.norm(temp_pca[:,0]),np.linalg.norm(temp_pca[:1])
        #meta info
        plt.title("K-means cluster {0:d}, mean={1:.3f}".format(i,mm))
        plt.xlabel("PC1, var={0:.3f}".format(temp_var[0]))
        plt.ylabel("PC2, var={0:.3f}".format(temp_var[1]))
        #transform data
        data_transformed = np.dot(data,temp_pca[:,:2]).transpose()
        #for i in xrange(data_transformed.shape[1]):
            #data_transformed[:,i]/=np.linalg.norm(data_transformed[:,i])
        #data_transformed = np.dot(temp_pca[:,:2].transpose(),data_transformed)
        #data_transformed = np.dot(data, temp_pca[:,:2])

        #data_transformed /= np.linalg.norm(data_transformed)
        print "after mean:", np.mean(data_transformed)
        #print data_transformed
        #print data_transformed
        #scatter plot data
        X=temp_pca[:,0].transpose()
        Y=temp_pca[:,1].transpose()
        #print title,i,zip(X,Y)
        plt.scatter(data_transformed[0,:],data_transformed[1,:])
        plt.quiver([0 for x in X], [0 for y in Y], X*mmax,Y*mmax, scale_units='xy', angles='xy', scale=1)
        for x,y,pfx in zip (X*mmax,Y*mmax,pfxs):
            plt.text(x+0.01,y+0.01, pfx)
        #quiver plot PCs
        #label the features
        #done
        #X+=temp_pca[0][:2]
        #Y+=temp_pca[1][:2]
        #C+=[colors[i],colors[i]]
        #plt.xscale('log')
        plt.xlim([-mmax,mmax])
        plt.ylim([-mmax,mmax])
        #plt.xlim([np.amin(X)*1.1,np.amax(X)*1.1])
        #plt.ylim([np.amin(Y)*1.1,np.amax(Y)*1.1])
        #plt.yscale('log')
    #plt.legend(["Cluster {0:d}".format(i) for i in xrange(K)])
    plt.savefig(title+" K_{0:d} ".format(K)+".png")

def plot (idxs, mean, pval, tags, pfxs, title):
    ch=chart.google_bubble(title = title,
            ylabel="negative log2 of p-value", xlabel="dNSAF mean (log scale)")
    DATA=[[] for i in pfxs]
    cols= ["Accession", "dNSAF","p-Value(log2)","K-means group"]
    for idx, m, p , tag  in  zip (idxs, mean, f(pval), tags):
        assert (len(m)==len(p)==len(pfxs))
        for i in xrange(len(pfxs)):
            if IDS[tag]!="N/A":
                DATA[i].append([IDS[tag], m[i], p[i],  str(idx)])
            else:
                DATA[i].append([tag, m[i], p[i],  str(idx)])
    print len(DATA[1])
    for pfx in pfxs:
        ch.load(cols, DATA[pfxs.index(pfx)], label=pfx+" Clustering Plot")
    ch.hlog="true"
    open(title+".html","w").write(ch.html())

def routine3():
    idxs = k(tags, baits_pvals, 4)
    SUBNETS["BAITS_PVAL_K4"]=(4,idxs)
    pca(4,idxs, baits_pvals,"Baits p-values",nat_pfxs[:-3])
    #pca(1,[0,0,0,0],[[1,2],[2,4],[-1,2],[-2,-4]],"test", ["X","Y=2X"])
    idxs = k(tags, f(baits_pvals), 4)
    SUBNETS["BAITS_PVAL_NL_K4"]=(4,idxs)
    #pca(4,idxs, f(baits_pvals),"Baits negaive log scale p-values",nat_pfxs[:-3])
    #plot(idxs,baits_means, baits_pvals, tags, nat_pfxs[:-3], "Baits_clusters_4")

    idxs = k(tags, grouped_pvals, 4)
    SUBNETS["GROUPED_PVAL_K4"]=(4,idxs)
    #pca(4,idxs, grouped_pvals,"Grouped p-values",nat_pfxs[-3:])
    idxs = k(tags, f(grouped_pvals), 4)
    SUBNETS["GROUPED_PVAL_NL_K4"]=(4,idxs)
    #pca(4,idxs, f(grouped_pvals),"Grouped negative log scale p-values",nat_pfxs[-3:])

    #plot(idxs,grouped_means, grouped_pvals, tags, nat_pfxs[-3:], "Grouped_clusters_4")
    #idxs = k(tags, baits_pvals, 3)
    #pca(3,idxs, baits_pvals,"Baits p-values",nat_pfxs[:-3])
    #idxs = k(tags, f(baits_pvals) , 3)
    #pca(3,idxs, f(baits_pvals),"Baits negative log scale p-values",nat_pfxs[:-3])

    #plot(idxs,baits_means, baits_pvals, tags, nat_pfxs[:-3], "Baits_clusters_3")
    #idxs = k(tags, grouped_pvals, 3)
    #pca(3,idxs, grouped_pvals,"Grouped p-values",nat_pfxs[-3:])
    #idxs = k(tags, f(grouped_pvals), 3)
    #pca(3,idxs, f(grouped_pvals),"Grouped negative log scale p-values",nat_pfxs[-3:])

    #plot(idxs,grouped_means, grouped_pvals, tags, nat_pfxs[-3:], "Grouped_clusters_3")
    #idxs = k(tags, baits_pvals, 5)
    #pca(5,idxs, baits_pvals,"Baits p-values",nat_pfxs[:-3])
    #idxs = k(tags, f(baits_pvals), 5)
    #pca(5,idxs, f(baits_pvals),"Baits negative log scale p-values",nat_pfxs[:-3])

    #plot(idxs,baits_means, baits_pvals, tags, nat_pfxs[:-3], "Baits_clusters_5")
    #idxs = k(tags, grouped_pvals, 5)
    #pca(5,idxs, grouped_pvals,"Grouped p-values",nat_pfxs[-3:])
    #idxs = k(tags,f( grouped_pvals), 5)
    #pca(5,idxs,f( grouped_pvals),"Grouped negative log scale p-values",nat_pfxs[-3:])
    #plot(idxs,grouped_means, grouped_pvals, ags, nat_pfxs[-3:], "Grouped_clusters_5")
    pkl.dump(SUBNETS,open("subnets.pkl","wb"))
#routine3()
def routine4():
    def gg(data, fname):
        ch=chart.jquery_datatable(title="Cluster Table")
        cols = ["Accession"]+["Cluster Index, K="+str(j) for j in range(20,51,5)]
        rows=[[i] for i in tags]
        for K in range(20, 51, 5):
            idxs=k(tags, data, K)
            rows=[i+[j] for i,j in zip(rows,idxs)]
        print len(cols), len(rows[0])
        ch.load(cols,rows)
        open("Kmeans_clusters_{0}.html".format(fname),"w").write(ch.html())

    gg ( baits_pvals, "baits_pvalues")
    gg ( grouped_pvals,  "groups_pvalues")
    gg ( baits_filtered_pvals,  "baits_filtered_pvalues")
    gg ( grouped_filtered_pvals,  "groups_filtered_pvalues")

def routine0():
    cols = ["Pt", "Mean","Exp", "Log(x+1)","x^(1/3)","x^0.5"]
    rows=[]
    ch=chart.google_scatter(title="Various Plots")

    rows = [[i,m,em,lm,m2,m5] for i,m,em,lm,m2,m5 in zip(xrange(len(MEANS)),sorted(MEANS.values()), sorted(map(math.exp,MEANS.values())),\
            sorted(map(lambda x: math.log(x+1,10),MEANS.values())), sorted(map(lambda x: x**(1/3), MEANS.values())),\
            sorted(map(lambda x: x**0.5, MEANS.values())))]
    ch.load(cols,rows)
    open("meanplots.html","w").write(ch.html())
routine0()

def routine5():
    def f1(col):
        return [[i**(0.5) for i in j] for j in col]
    def f2(col):
        return [[math.exp(i) for i in j] for j in col]
    KK=[20,35,50]
    ch=chart.jquery_datatable(title="Cluster Table")
    cols = ["Accession","Name"]+["Cluster Index, K={0:d}, data = {1}".format(j,i)\
            for j in KK\
            for i in ["raw pval","filtered_pvals","negative log10 of p values","simple mean","filtered_mean",\
            "simple mean, square root","filtered mean, square root","simple mean, exp", "filtered mean, exp"]\
            ]
    rows =[[i,j] for i,j in zip(tags,htags)]
    for K in KK:
        P=k(tags,baits_pvals, K)
        FP=k(tags,baits_filtered_pvals, K)
        NP=k(tags,f(baits_pvals), K)
        M=k(tags,baits_means, K)
        FM=k(tags,baits_filtered_means, K)
        F1M= k(tags,f1(baits_means), K)
        F1FM= k(tags,f1(baits_filtered_means), K)
        F2M= k(tags,f2(baits_means), K)
        F2FM= k(tags,f2(baits_filtered_means), K)
        rows=[rr+ [p, fp,np, m, fm , f1m, f1fm, f2m, f2fm] for rr, p, fp,np, m, fm ,f1m, f1fm, f2m, f2fm in zip(rows, P,FP,NP,M,FM,F1M, F1FM, F2M, F2FM)]
    ch.load(cols,rows)
    open("clusters_various.html","w").write(ch.html())
    #get multiple data and pfxs and generate the k means clusters table

def routine5():
    def f1(col):
        return [[i**(0.5) for i in j] for j in col]
    def f2(col):
        return [[math.exp(i) for i in j] for j in col]
    KK=[20,35,50]
    ch=chart.jquery_datatable(title="Cluster Table")
    cols = ["Accession","Name"]+["Cluster Index, K={0:d}, data = {1}".format(j,i)\
            for j in KK\
            for i in ["raw pval","filtered_pvals","negative log10 of p values","simple mean","filtered_mean",\
            "simple mean, square root","filtered mean, square root","simple mean, exp", "filtered mean, exp"]\
            ]
    rows =[[i,j] for i,j in zip(tags,htags)]
    for K in KK:
        P=k(tags,baits_pvals, K)
        FP=k(tags,baits_filtered_pvals, K)
        NP=k(tags,f(baits_pvals), K)
        M=k(tags,baits_means, K)
        FM=k(tags,baits_filtered_means, K)
        F1M= k(tags,f1(baits_means), K)
        F1FM= k(tags,f1(baits_filtered_means), K)
        F2M= k(tags,f2(baits_means), K)
        F2FM= k(tags,f2(baits_filtered_means), K)
        rows=[rr+ [p, fp,np, m, fm , f1m, f1fm, f2m, f2fm] for rr, p, fp,np, m, fm ,f1m, f1fm, f2m, f2fm in zip(rows, P,FP,NP,M,FM,F1M, F1FM, F2M, F2FM)]
    ch.load(cols,rows)
    open("clusters_various.html","w").write(ch.html())
    #get multiple data and pfxs and generate the k means clusters table
def routine6():
    def f1(col):
        return [[i**(0.5) for i in j] for j in col]
    def f2(col):
        return [[math.exp(i) for i in j] for j in col]
    KK=range(20,80,5)
    cols = ["Accession","Name"]+["Cluster Index, K={0:d}, data = {1}".format(j,i)\
            for i in ["raw pval","filtered_pvals"]\
            for j in KK ]
    print cols
    d = tabular()
    d<< (tags, "Accessions")
    d<< (htags, "Names")
    ii = 2
    for _d in [baits_pvals, baits_filtered_pvals]:
        for K in KK:
            print len(d.data),d.shape(),len(_d)
            d<< (k(tags, _d, K, tries=200),cols[ii])
            ii+=1

    print len (d.data), len(d.data[0])
    open("pvalclusters.xls","w").write(repr(d))
routine6()

