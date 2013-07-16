#include<iostream>
#include<boost/python.hpp>
#include<ap.h>

void hello()
{
    std::cout<<"Hi There, this is a boost/python test program.\n";
    std::cout<<"Now testing python extension.\n";
    return;
}
BOOST_PYTHON_MODULE(test_ext)
{
    using namespace boost::python;
    def("hello", hello);
}
//int main(int argc, const char** argv)
//{
//    hello();
//    return 0;
//}
