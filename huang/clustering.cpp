#include<iostream>
#include<boost/python.hpp>
#include"alglib/ap.h"
#include"alglib/dataanalysis.h"

boost::python::list svdpca (boost::python::str ss,
        int nsams, int nvars)
// Wrapper for alglib pca
// assuming ss is formatted correctly as a 2d array
{
    const char* s = boost::python::extract<const char*> (ss);
    alglib::real_2d_array arr(s);
    //long nnsams = boost::python::extract<long> (nsams);
    //long nnvars = boost::python::extract<long> (nvars);
    //type conversion
    alglib::ae_int_t info=0;
    alglib::real_1d_array Var;
    alglib::real_2d_array Bas;
    alglib::pcabuildbasis (s, nsams, nvars, info, Var, Bas);
    if (info==0) std::cerr<<"Unknown Error\n";
    else if (info==-4) std::cerr<<"No conversion\n";
    else if (info==-1) std::cerr<<"Wrong format\n";
    else std::cerr<<"SVD PCA success\n";
    boost::python::str out1=boost::python::str(Var.tostring(3));
    boost::python::str out2=boost::python::str(Bas.tostring(3));
    boost::python::list out;
    out.append(out1);
    out.append(out2);
    return out;
}

boost::python::str kmeans (boost::python::str ss,
        int K, int resets = 5)
{
    using namespace alglib;
    const char* sss = boost::python::extract<const char*> (ss);
    clusterizerstate s;
    kmeansreport rep;
    real_2d_array arr(sss) ;
    clusterizersetpoints(s, arr, 2);
    //must be euclidean distances
    clusterizersetkmeanslimits(s, resets, 0);
    //set a low random reset number of 5 times
    clusterizerrunkmeans(s, K, rep);
    int flag = int(rep.terminationtype);
    if (flag == 1) std::cerr<<"KMeans Clustering Successful!\n";
    else if (flag==-4) std::cerr<<"No convergence\n";
    else if (flag==-1) std::cerr<<"Wrong format\n";
    else std::cerr<<"Unknown Error type: "<<flag<<std::endl;
    return boost::python::str(rep.cidx.tostring());
}
BOOST_PYTHON_MODULE(clustering)
{
    using namespace boost::python;
    def ("svdpca",svdpca);
    def ("kmeans",kmeans);
}

