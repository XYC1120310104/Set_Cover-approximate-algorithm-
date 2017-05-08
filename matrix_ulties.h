#ifndef __MATRIX_ULTIES_H__
#include <vector>
using namespace std;
vector<vector<double> > matrix_multi(vector<vector<double> > &a,vector<vector<double> > &b,int m,int n,int q);
void matrix_num_multi(vector<vector<double> > &a,double number,int m,int n);
vector<vector<double> > matrix_sub(vector<vector<double> > &a,vector<vector<double> > &b,int m,int n);
vector<vector<double> > matrix_div(vector<vector<double> > &a,vector<vector<double> > &b,int m,int n);
vector<vector<double> > matrix_plus(vector<vector<double> > &a,vector<vector<double> > &b,int m,int n);
vector<vector<double> > diag(vector<vector<double> > &a,int m,int n=1);//n==1
vector<vector<double> > matrix_trans(vector<vector<double> > &a,int m,int n);

vector<vector<double> > matrix_multi(vector<vector<double> > &a,vector<vector<double> > &b,int m,int n,int q)//a*b
{
    vector<vector<double> > c(m,vector<double>(q));
    int i,j,k;
    for(i=0;i<m;i++)
        for(j=0;j<q;j++)
        {
            c[i][j]=0;
            for(k=0;k<n;k++)
            c[i][j]=c[i][j]+a[i][k]*b[k][j];
        }
    return c;
}

void matrix_num_multi(vector<vector<double> > &a,double number,int m,int n)//number*a
{
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            a[i][j]=number*a[i][j];
        }
    }
}

vector<vector<double> > matrix_sub(vector<vector<double> > &a,vector<vector<double> > &b,int m,int n)//a-b
{
    vector<vector<double> > c(m,vector<double>(n));
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
        {
            c[i][j]=a[i][j]-b[i][j];
        }
    return c;
}

vector<vector<double> > matrix_div(vector<vector<double> > &a,vector<vector<double> > &b,int m,int n)//a./b
{
    vector<vector<double> > res(m,vector<double>(n));
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            res[i][j]=a[i][j]/b[i][j];
        }
    }
    return res;
}

vector<vector<double> > matrix_plus(vector<vector<double> > &a,vector<vector<double> > &b,int m,int n)//a+b
{
    vector<vector<double> > c(m,vector<double>(n));
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
        {
            c[i][j]=a[i][j]+b[i][j];
        }
    return c;
}

vector<vector<double> > diag(vector<vector<double> > &a,int m,int n)
{
    vector<vector<double> > res(m,vector<double>(m));
    for (int i=0;i<m;i++){
        res[i][i]=a[i][0];
    }
    return res;
}

vector<vector<double> > matrix_trans(vector<vector<double> > &a,int m,int n)
{
    vector<vector<double> > res(n,vector<double>(m,0));
    for (int i=0;i<n;i++){
        for (int j=0;j<m;j++){
            res[i][j]=a[j][i];
        }
    }
    return res;
}
#endif // __MATRIX_ULTIES_H__
