import os, sys, string, numpy, math
from copy import deepcopy
import cPickle as pkl
#import crick.viz.js.google_chart_js as chart
import crick.annotate.cybert.crick_cybert_annotator as anno
import crick.builders.protein_protein_edges.StringBuilder as sb
from crick.utils.tabular import tabular
from crick.analyzer import ProteinAnalyzer
from crick.const import species
from wrapper import convertIDs
targets = os.listdir("./csv")
print targets
tags = open("./csv/mean_all.txt","r").read().split("\n")
tags = sorted([i.strip() for i in tags])
print tags
ids = convertIDs (fromType="UNIPROT_ACCESSION", toType="OFFICIAL_GENE_SYMBOL",ids=list(tags))
IDS = dict(ids)
print len(tags)
htags= [IDS[i] if IDS[i]!="N/A" else i for i in tags]
IDS = dict(zip(tags,htags))
PVALS=pkl.load(open("./pkl/NatDatapvals_all_all.pkl","rb"))
MEANS=pkl.load(open("./pkl/NatDatameans_all_all.pkl","rb"))

BAITS=["Pre1","Pre10","Rpt1","Rpn1","Rpt5","Rpt6","Rpn11"]+["Base","Lid","20S"]
nat_pfxs= ["Nat_"+i for i in BAITS]
den_pfxs= ["Den_"+i for i in BAITS]
#nat_gpfxs= ["Nat_"+i for i in GBAITS]
#den_gpfxs= ["Den_"+i for i in GBAITS]
SUBNETS=pkl.load(open("subnets.pkl","rb"))
#from here everything is native
all_means = [[float(MEANS[(tag,pfx)]) if (tag,pfx) in MEANS else 0.0 for pfx in nat_pfxs] for tag in tags]
all_pvals = [[float(PVALS[(tag,pfx)]) if (tag,pfx) in PVALS else 1.0 for pfx in nat_pfxs] for tag in tags]
#print all_pvals
#print all_pvals
sig_flags = [[float(i)<0.01 for i in j] for j in all_pvals]
sig_flags = [[float(i)<0.05 for i in j] for j in all_pvals]
#print sig_flags

filtered_means=[]
filtered_pvals=[]
#now produce filtered data

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
baits_pvals = [i[:-3] for i in all_pvals]
baits_filtered_means = [i[:-3] for i in filtered_means]
baits_filtered_pvals = [i[:-3] for i in filtered_pvals]
print nat_pfxs[:-3], nat_pfxs[-3:]

grouped_means = [i[-3:] for i in all_means]
grouped_pvals = [i[-3:] for i in all_pvals]
grouped_filtered_means = [i[-3:] for i in filtered_means]
grouped_filtered_pvals = [i[-3:] for i in filtered_pvals]
def routine2():
    #build networks
    an = anno.crick_tester(species.hg19)
    an.n.add_proteins(tags)
    an.name="huang_test"
    an.closed_ppi()
    an.do_annotation_routine()
    an.exportfig("huang_grid_closed")
    an.pickle("huang_grid")
    del an
    an = anno.crick_tester(species.hg19)
    an.n.add_proteins(tags)
    an.name="huang_test"
#string_itype="String_Functional_Interaction"
    sbuilder=sb.StringProteinProteinEdgeBuilder(an.n)
    sbuilder.build_network(closed_network=True,string_itype="binding", string_score=600);
    an.do_annotation_routine()
    an.exportfig("huang_string_closed_all")
    an.pickle("huang_string")
    del an
    an = anno.crick_tester(species.hg19)
    an.n.add_proteins(tags)
    an.name="huang_test"
#string_itype="String_Functional_Interaction"
    sbuilder=sb.StringProteinProteinEdgeBuilder(an.n)
    sbuilder.build_network(closed_network=True,string_itype="binding", string_score=600);
    an.do_annotation_routine()
    an.exportfig("huang_string_closed_600")
    an.pickle("huang_string_600")
#routine2()
def routine4 ():
    for fname, chunk in SUBNETS.items():
        for i in xrange (chunk[0]):
            tmp_tags=  filter(lambda x: chunk[1][tags.index(x)]==i,tags)
            an = anno.crick_tester(species.hg19)
            an.n.add_proteins(tmp_tags)
            #an.closed_ppi()
            #an.do_annotation_routine(f=tmp_tags+["Kmeans-cluster{0:d}".format(i)])
            #an.exportfig(fname+"_GRID_{0:d}_of_{1:d}".format(i,chunk[0]))
            #an.pickle(fname+"_GRID_{0:d}_of_{1:d}".format(i,chunk[0]))
            sbuilder=sb.StringProteinProteinEdgeBuilder(an.n)
            sbuilder.build_network(closed_network=True,string_itype="binding", string_score=600);
            an.do_annotation_routine(f=tmp_tags+["Kmeans-cluster{0:d}".format(i)])
            an.exportfig(fname+"_STRING_600_{0:d}_of_{1:d}".format(i,chunk[0]))
            an.pickle(fname+"_STRING_600_{0:d}_of_{1:d}".format(i,chunk[0]))
def routine5():
    an=anno.crick_tester()
    an=an.unpickle("/home/yul13/tmp/initial.pkl")
    ap = ProteinAnalyzer(an.n,"./analysis")
    #ap.DAVIDAnalysis()
    tt= tabular()
    tt.loadf("./analysis/DAVIDReports/EnrichmentReport.xls")
    tt.col_map(lambda x:",".join([IDS[i] for i in x.split(", ")]),"geneIds")
    open("Enrichment.csv","w").write(repr(tt))

def routine6():
    K=30
    #idxs = t % "Cluster Index, K=30, data = filtered_pvals"
    #print idxs

    t= tabular()
    t.loadf("pvalclusters.xls")
    t.order()
    print "O list is now" , t.order_list
    t.col_map(int, "Cluster Index, K=30, data = filtered_pvals")
    for ii in xrange(K):
        print "O list is now" , t.order_list
        ns = t.filter(lambda x: int(x)==ii,"Cluster Index, K=30, data = filtered_pvals" )
        print t.order_list , len(ns.data), (ns % 0), len(t.data), len (t.data[0])
        print ns % "Cluster Index, K=30, data = filtered_pvals"
        ns= (ns % 0)
        print "got {0:d} proteins for i={1:d}".format(len(ns),ii)
        t.reorder()
        try:
            an = anno.crick_tester(species.hg19)
            an.n.add_proteins(ns)
            an.closed_ppi()
            an.do_annotation_routine(f=ns+["Kmeans Cluster{0:d} of 30".format(ii)])
            an.exportfig("Huang_K30_{0:d}.xgmml".format(ii))
            del an

            an = anno.crick_tester(species.hg19)
            an.n.add_proteins(ns)
            sbuilder=sb.StringProteinProteinEdgeBuilder(an.n)
            sbuilder.build_network(closed_network=True,string_itype="binding", string_score=600);
            an.do_annotation_routine(f=ns+["Kmeans-cluster{0:d} of 30".format(ii)])
            an.exportfig("Huang_STRING_600_K30_{0:d}".format(ii))

            del an
            del sbuilder
            an = anno.crick_tester(species.hg19)
            an.n.add_proteins(ns)
            sbuilder=sb.StringProteinProteinEdgeBuilder(an.n)
            sbuilder.build_network(closed_network=True,string_itype="binding", string_score=900);
            an.do_annotation_routine(f=ns+["Kmeans-cluster{0:d} of 30".format(ii)])
            an.exportfig("Huang_STRING_900_K30_{0:d}".format(ii))

            ap = ProteinAnalyzer(an.n,"./tmpanalysis")
            ap.DAVIDAnalysis()
            tt= tabular()
            tt.loadf("./tmpanalysis/DAVIDReports/EnrichmentReport.xls")
            tt.col_map(lambda x:",".join([IDS[i] for i in x.split(", ")]),"geneIds")
        except:
            continue
        open("subnetworks/analysis/Huang_K30_{0:d}.csv".format(ii),"w").write(repr(tt))
#routine6()
def routine7():
    t= tabular()
    t.loadf("pvalclusters.xls")
    t.order()
    print "O list is now" , t.order_list
    t.col_map(int, "Cluster Index, K=30, data = filtered_pvals")

    r = tabular()
    r.loadf("./allComplexes.csv", sep=";")

    prots =  r % "subunits (UniProt IDs)"
    CPX = {}
    K=30
    for i in prots:
        for j in i.split(","):
            j=j.strip("(").strip(")")
            CPX[j]=prots.index(i)
    for ii in xrange(K):
        ns = t.filter(lambda x: int(x)==ii,"Cluster Index, K=30, data = filtered_pvals" )
        #print t.order_list , len(ns.data), (ns % 0), len(t.data), len (t.data[0])
        #print ns % "Cluster Index, K=30, data = filtered_pvals"
        ns= (ns % 0)
        print "got {0:d} proteins for i={1:d}".format(len(ns),ii)
        rr = tabular()
        rr.header =deepcopy (r.header)
        rr.nc = r.nc
        rr << ([],"Number of Subunits in cluster")
        print rr
        hits = {}
        for i in ns:
            if i in CPX:
                print "Found a complex:", (r % 1)[CPX[i]]
                if CPX[i] in hits:
                    hits [CPX[i]]+=1
                else:
                    hits [CPX[i]]=1
        print hits
        for idx, num_hit in hits.items():
            rr < (r[idx]+[num_hit])
        rr.order(f=lambda x: int(x))
        open("./subnetworks/complex/Cluster_{0:d}.xls".format(ii),"w").write(repr(rr))
        t.reorder()
routine7()
