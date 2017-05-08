#ifndef __MATRIX_INV_H__
#define __MATRIX_INV_H__

#include <vector>
using namespace std;
void temp(vector<double> &aa,vector<double> &bb,int n);
double fun(vector<vector<double> > &array,int n);//求矩阵行列式
double matrix_inv(vector<vector<double> > &a,vector<vector<double> > &c,int n);//矩阵求逆

void temp(vector<double> &aa,vector<double> &bb,int n)
{
    for(int i=0; i<n; i++){
        double temp1=aa[i];
        aa[i]=bb[i];
        bb[i]=temp1;
    }
}

double fun(vector<vector<double> > &array,int n)
{
    int ii,jj,k,u;
    double det1=1,yin;
    for(ii=0; ii<n; ii++)
    {
        if(array[ii][ii]==0)
            for(jj=ii; jj<n; jj++)
            {
                if(array[jj][ii]!=0)
                    temp(array[ii],array[jj],n);
            }
        for(k=ii+1; k<n; k++)
        {
            yin=-1*array[k][ii]/array[ii][ii];
            for(u=0; u<n; u++)
            {
                array[k][u]=array[k][u]+array[ii][u]*yin;
            }
        }
    }
    for(ii=0; ii<n; ii++)
        det1=det1*array[ii][ii];
    return (det1);
}
double matrix_inv(vector<vector<double> > &a,vector<vector<double> > &c,int n)//c为a逆矩阵
{
    vector<vector<double> > b(n,vector<double>(2*n));
    double det1,yinzhi;
    double bb;
    int i,j,k,u;
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            b[i][j]=a[i][j];
    for(j=0; j<n; j++)
        b[j][n+j]=1;
    det1=fun(a,n);
    for(i=0; i<n; i++)
    {
        if(b[i][i]==0)
            for(j=i; j<n; j++)
            {
                if(b[j][i]!=0) temp(b[i],b[j],n);
            }
        for(k=i+1; k<n; k++)
        {
            yinzhi=-1*b[k][i]/b[i][i];
            for(u=0; u<2*n; u++)
            {
                b[k][u]=b[k][u]+b[i][u]*yinzhi;
            }
        }
    }
    det1=fun(a,n);
    if(det1==0)      return det1;
    if(det1!=0)
    {
        for(i=0; i<n; i++)
        {
            bb=b[i][i];
            for(j=0; j<2*n; j++)
                b[i][j]=b[i][j]/bb;
        }
        for(i=n-1; i>0; i--)
            for(k=0; k<i; k++)
            {
                bb=b[k][i];
                for(u=0; u<2*n; u++)
                    b[k][u]=b[k][u]-bb*b[i][u];
            }
    }
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            c[i][j]=b[i][j+n];
    return det1;
}
#endif // __MATRIX_INV__H__



