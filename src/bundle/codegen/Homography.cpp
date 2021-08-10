#include <cmath>
double beta(double* c, double* b, double* d, double* p1, double* p2, double* p3, double* p4, double* a)
{
 return (c[0]*d[3]*std::sin(a[0])*std::cos(a[1]) - c[0]*d[3]*std::sin(a[1])*std::cos(a[0]) - c[0]*d[6]*u2,1*std::sin(a[0])*std::sin(a[1]) + c[0]*d[6]*u2,2*std::sin(a[0])*std::sin(a[1]) + c[0]*d[6]*v2,1*std::sin(a[1])*std::cos(a[0]) - c[0]*d[6]*v2,2*std::sin(a[0])*std::cos(a[1]) - c[3]*d[0]*std::sin(a[0])*std::cos(a[1]) + c[3]*d[0]*std::sin(a[1])*std::cos(a[0]) + c[3]*d[6]*u2,1*std::sin(a[0])*std::cos(a[1]) - c[3]*d[6]*u2,2*std::sin(a[1])*std::cos(a[0]) - c[3]*d[6]*v2,1*std::cos(a[0])*std::cos(a[1]) + c[3]*d[6]*v2,2*std::cos(a[0])*std::cos(a[1]) + c[6]*d[0]*u2,1*std::sin(a[0])*std::sin(a[1]) - c[6]*d[0]*u2,2*std::sin(a[0])*std::sin(a[1]) - c[6]*d[0]*v2,1*std::sin(a[1])*std::cos(a[0]) + c[6]*d[0]*v2,2*std::sin(a[0])*std::cos(a[1]) - c[6]*d[3]*u2,1*std::sin(a[0])*std::cos(a[1]) + c[6]*d[3]*u2,2*std::sin(a[1])*std::cos(a[0]) + c[6]*d[3]*v2,1*std::cos(a[0])*std::cos(a[1]) - c[6]*d[3]*v2,2*std::cos(a[0])*std::cos(a[1]))/(b[0]*c[3]*std::sin(a[0])*std::cos(a[1]) - b[0]*c[3]*std::sin(a[1])*std::cos(a[0]) - b[0]*c[6]*u2,1*std::sin(a[0])*std::sin(a[1]) + b[0]*c[6]*u2,2*std::sin(a[0])*std::sin(a[1]) + b[0]*c[6]*v2,1*std::sin(a[1])*std::cos(a[0]) + b[0]*c[6]*v2,2*std::sin(a[0])*std::cos(a[1]) - b[3]*c[0]*std::sin(a[0])*std::cos(a[1]) - b[3]*c[0]*std::sin(a[1])*std::cos(a[0]) + b[3]*c[6]*u2,1*std::sin(a[0])*std::cos(a[1]) - b[3]*c[6]*u2,2*std::sin(a[1])*std::cos(a[0]) - b[3]*c[6]*v2,1*std::cos(a[0])*std::cos(a[1]) + b[3]*c[6]*v2,2*std::cos(a[0])*std::cos(a[1]) + b[6]*c[0]*u2,1*std::sin(a[0])*std::sin(a[1]) - b[6]*c[0]*u2,2*std::sin(a[0])*std::sin(a[1]) - b[6]*c[0]*v2,1*std::sin(a[1])*std::cos(a[0]) + b[6]*c[0]*v2,2*std::sin(a[0])*std::cos(a[1]) - b[6]*c[3]*u2,1*std::sin(a[0])*std::cos(a[1]) + b[6]*c[3]*u2,2*std::sin(a[1])*std::cos(a[0]) + b[6]*c[3]*v2,1*std::cos(a[0])*std::cos(a[1]) - b[6]*c[3]*v2,2*std::cos(a[0])*std::cos(a[1]));
}

double gamma(double* c, double* b, double* d, double* p1, double* p2, double* p3, double* p4, double* a)
{
 return (-b[0]*d[3]*std::sin(a[0])*std::cos(a[1]) + b[0]*d[3]*std::sin(a[1])*std::cos(a[0]) + b[0]*d[6]*u2,1*std::sin(a[0])*std::sin(a[1]) - b[0]*d[6]*u2,2*std::sin(a[0])*std::sin(a[1]) - b[0]*d[6]*v2,1*std::sin(a[1])*std::cos(a[0]) + b[0]*d[6]*v2,2*std::sin(a[0])*std::cos(a[1]) - b[3]*d[0]*std::sin(a[0])*std::cos(a[1]) - b[3]*d[0]*std::sin(a[1])*std::cos(a[0]) - b[3]*d[6]*u2,1*std::sin(a[0])*std::cos(a[1]) + b[3]*d[6]*v2,1*std::cos(a[0])*std::cos(a[1]) + b[3]*d[6]*v2,2*std::sin(a[1])*std::cos(a[0]) + b[3]*d[6]*v2,2*std::cos(a[0])*std::cos(a[1]) - b[6]*d[0]*u2,1*std::sin(a[0])*std::sin(a[1]) + b[6]*d[0]*u2,2*std::sin(a[0])*std::sin(a[1]) - b[6]*d[0]*v2,1*std::sin(a[1])*std::cos(a[0]) - b[6]*d[0]*v2,2*std::sin(a[0])*std::cos(a[1]) + b[6]*d[3]*u2,1*std::sin(a[0])*std::cos(a[1]) - b[6]*d[3]*u2,2*std::sin(a[1])*std::cos(a[0]) - b[6]*d[3]*v2,1*std::cos(a[0])*std::cos(a[1]) + b[6]*d[3]*v2,2*std::cos(a[0])*std::cos(a[1]))/(b[0]*c[3]*std::sin(a[0])*std::cos(a[1]) - b[0]*c[3]*std::sin(a[1])*std::cos(a[0]) - b[0]*c[6]*u2,1*std::sin(a[0])*std::sin(a[1]) + b[0]*c[6]*u2,2*std::sin(a[0])*std::sin(a[1]) + b[0]*c[6]*v2,1*std::sin(a[1])*std::cos(a[0]) - b[0]*c[6]*v2,2*std::sin(a[0])*std::cos(a[1]) - b[3]*c[0]*std::sin(a[0])*std::cos(a[1]) + b[3]*c[0]*std::sin(a[1])*std::cos(a[0]) + b[3]*c[6]*u2,1*std::sin(a[0])*std::cos(a[1]) - b[3]*c[6]*u2,2*std::sin(a[1])*std::cos(a[0]) - b[3]*c[6]*v2,1*std::cos(a[0])*std::cos(a[1]) + b[3]*c[6]*v2,2*std::cos(a[0])*std::cos(a[1]) + b[6]*c[0]*u2,1*std::sin(a[0])*std::sin(a[1]) - b[6]*c[0]*u2,2*std::sin(a[0])*std::sin(a[1]) - b[6]*c[0]*v2,1*std::sin(a[1])*std::cos(a[0]) + b[6]*c[0]*v2,2*std::sin(a[0])*std::cos(a[1]) - b[6]*c[3]*u2,1*std::sin(a[0])*std::cos(a[1]) + b[6]*c[3]*u2,2*std::sin(a[1])*std::cos(a[0]) + b[6]*c[3]*v2,1*std::cos(a[0])*std::cos(a[1]) - b[6]*c[3]*v2,2*std::cos(a[0])*std::cos(a[1]));
}

double nullSpaceVecB0(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(-p1[0]*p2[0]*p3[1] + p1[0]*p2[0]*p5[1] + p1[1]*p3[0]*p4[0] - p1[1]*p5[0]*p6[0] - p3[0]*p4[0]*p5[1] + p3[1]*p5[0]*p6[0])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecB1(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(p1[0]*p2[0]*p3[0] - p1[0]*p2[0]*p5[0] - p1[0]*p3[0]*p4[0] + p1[0]*p5[0]*p6[0] + p3[0]*p4[0]*p5[0] - p3[0]*p5[0]*p6[0])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecB2(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(-p1[0]*p2[0]*p3[0]*p5[1] + p1[0]*p2[0]*p3[1]*p5[0] + p1[0]*p3[0]*p4[0]*p5[1] - p1[0]*p3[1]*p5[0]*p6[0] - p1[1]*p3[0]*p4[0]*p5[0] + p1[1]*p3[0]*p5[0]*p6[0])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecB3(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(-p1[0]*p2[1]*p3[1] + p1[0]*p2[1]*p5[1] + p1[1]*p3[0]*p4[1] - p1[1]*p5[0]*p6[1] - p3[0]*p4[1]*p5[1] + p3[1]*p5[0]*p6[1])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecB4(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(p1[0]*p2[1]*p3[0] - p1[0]*p2[1]*p5[0] - p1[0]*p3[0]*p4[1] + p1[0]*p5[0]*p6[1] + p3[0]*p4[1]*p5[0] - p3[0]*p5[0]*p6[1])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecB5(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(-p1[0]*p2[1]*p3[0]*p5[1] + p1[0]*p2[1]*p3[1]*p5[0] + p1[0]*p3[0]*p4[1]*p5[1] - p1[0]*p3[1]*p5[0]*p6[1] - p1[1]*p3[0]*p4[1]*p5[0] + p1[1]*p3[0]*p5[0]*p6[1])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecB6(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 1;
}

double nullSpaceVecB7(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecB8(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecC0(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(-p1[1]*p2[0]*p3[1] + p1[1]*p2[0]*p5[1] + p1[1]*p3[1]*p4[0] - p1[1]*p5[1]*p6[0] - p3[1]*p4[0]*p5[1] + p3[1]*p5[1]*p6[0])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecC1(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(-p1[0]*p3[1]*p4[0] + p1[0]*p5[1]*p6[0] + p1[1]*p2[0]*p3[0] - p1[1]*p2[0]*p5[0] - p3[0]*p5[1]*p6[0] + p3[1]*p4[0]*p5[0])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecC2(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(p1[0]*p3[1]*p4[0]*p5[1] - p1[0]*p3[1]*p5[1]*p6[0] - p1[1]*p2[0]*p3[0]*p5[1] + p1[1]*p2[0]*p3[1]*p5[0] + p1[1]*p3[0]*p5[1]*p6[0] - p1[1]*p3[1]*p4[0]*p5[0])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecC3(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(-p1[1]*p2[1]*p3[1] + p1[1]*p2[1]*p5[1] + p1[1]*p3[1]*p4[1] - p1[1]*p5[1]*p6[1] - p3[1]*p4[1]*p5[1] + p3[1]*p5[1]*p6[1])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecC4(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(-p1[0]*p3[1]*p4[1] + p1[0]*p5[1]*p6[1] + p1[1]*p2[1]*p3[0] - p1[1]*p2[1]*p5[0] - p3[0]*p5[1]*p6[1] + p3[1]*p4[1]*p5[0])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecC5(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return -(p1[0]*p3[1]*p4[1]*p5[1] - p1[0]*p3[1]*p5[1]*p6[1] - p1[1]*p2[1]*p3[0]*p5[1] + p1[1]*p2[1]*p3[1]*p5[0] + p1[1]*p3[0]*p5[1]*p6[1] - p1[1]*p3[1]*p4[1]*p5[0])/(p1[0]*p3[1] - p1[0]*p5[1] - p1[1]*p3[0] + p1[1]*p5[0] + p3[0]*p5[1] - p3[1]*p5[0]);
}

double nullSpaceVecC6(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecC7(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 1;
}

double nullSpaceVecC8(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecD0(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecD1(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecD2(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecD3(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecD4(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecD5(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecD6(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecD7(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 0;
}

double nullSpaceVecD8(double* p1, double* p2, double* p3, double* p4, double* p5, double* p6)
{
 return 1;
}

