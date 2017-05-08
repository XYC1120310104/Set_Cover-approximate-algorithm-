#include <iostream>
#include <set>
#include <vector>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <fstream>

//#include "matrix_inverse.h"
//#include "matrix_ulties.h"

using namespace std;

void print_F(vector<set<int> > F,string file);
void print_S(set<int> S,string file);
void print_matrix(vector<vector<double> > &A,string file);
int M=20;
set<int> X_generalization(int N){
    set<int> X;
    for (int i=0;i<N;i++){
        X.insert(i);
    }
    return X;
}
void removeDup(set<int> &S,set<int> &X){
    set<int>::iterator iter;
    for (iter=S.begin();iter!=S.end();iter++){
        if (X.count(*iter)) X.erase(*iter);
    }
}
set<int> random_select_element(int M,set<int> &X,bool flag=true){
    set<int>::iterator iter;
    set<int> S;
    int N=X.size(),i=0;
    for (iter=X.begin();iter!=X.end();iter++){
        double p_select=M/(double)(N-i),p_rand=rand()/(double)(RAND_MAX);
        if (p_select>=p_rand){
            S.insert(*iter);
            M--;
        }
        i++;
        if (M==0) break;
    }
    if (flag) removeDup(S,X);
    return S;
}

set<int> addAll(set<int> &U_S,set<int> S){
    set<int>::iterator S_iter=S.begin();
    for (S_iter=S.begin();S_iter!=S.end();S_iter++){
        U_S.insert(*S_iter);
    }
    return U_S;
}

vector<set<int> > F_generalization(int N,set<int> &X){
    vector<set<int> > F;
    F.push_back(random_select_element(M,X,true));
    int i=0;
    set<int> U_S;
    while(X.size()>=20){
        int n=2+(M-1)*(rand()/(double)(RAND_MAX));//生成n,2到M之间
        int x=1+(n-1)*(rand()/(double)(RAND_MAX));//生成x,1到n之间
        set<int> S=random_select_element(x,X,true);
        U_S=addAll(U_S,F[i]);
        S=addAll(S,random_select_element(n-x,U_S,true));
        F.push_back(S);
        i++;
    }
    F.push_back(X);
    //print_S(U_S);
    int F_num=F.size();
    for (int i=0;i<(N-F_num);i++){
        int n=2+(M-2)*(rand()/(double)(RAND_MAX));
        set<int> S=random_select_element(n,X,false);
        F.push_back(S);
    }
    return F;
}

int intersection(set<int> &S,set<int> &X){//求交运算
    set<int>::iterator iter;
    int cnt=0;
    for (iter=S.begin();iter!=S.end();iter++){
        if (X.count(*iter)) cnt++;
    }
    return cnt;
}

int set_sub(set<int> &S,set<int> &X,vector<double> &price){
    set<int>::iterator iter;
    int cnt=0;
    for (iter=S.begin();iter!=S.end();iter++){
        if (X.count(*iter)) cnt++;
    }
    cnt=S.size()-cnt;
    for (iter=S.begin();iter!=S.end();iter++){
        if (!X.count(*iter)){
            price[*iter]=1.0/(double)cnt;
        }
    }
    return cnt;
}

set<int> max_element_cover(vector<set<int> > &F,set<int> &X){
    int num=F.size();
    set<int> res;
    int max_cover=INT_MIN,index=0;
    for (int i=0;i<num;i++){
        int insec=intersection(F[i],X);
        if (max_cover<insec) {
            res=F[i];
            max_cover=insec;
            index=i;
        }
    }
    F.erase(F.begin()+index);
    removeDup(res,X);
    return res;
}

vector<set<int> > Greedy_Set_Cover(vector<set<int> > F,set<int> X){
    vector<set<int> > res;
    while(!X.empty()){
        set<int> S=max_element_cover(F,X);
        res.push_back(S);
    }
    return res;
}
//next solution is LP:affine-scaling
// min c'x
//    Ax>=1
//    Ex>=0
//    Ex<=1
//normalization:
//min c'x
//  -Ax<=-1
//  -Ex<=0
//   Ex<=1
/*
double get_positive_min_value(vector<vector<double> > &vk,vector<vector<double> > &hv,int m,int n){
    double res=INT_MAX;
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            if (hv[i][j]<0){
                double tmp=(vk[i][j]/(-1.0*hv[i][j]));
                res=res>tmp?tmp:res;
            }
        }
    }
    return res;
}

double get_min_value(vector<vector<double> > &matrix,int m,int n){
    double res=INT_MAX;
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            res=res>matrix[i][j]?matrix[i][j]:res;
        }
    }
    return res;
}
int get_LP_parameters(vector<set<int> > &F,set<int> &X,vector<vector<double> > &A,vector<vector<double> > &b,int N){
    set<int>::iterator iter;
    int i=0,f=0;
    for (iter=X.begin();iter!=X.end();iter++){
        int tmp=0;
        for (int j=0;j<N;j++){
            if (F[j].count(*iter)){
                A[j][i]=-1;
                tmp++;
            }
        }
        f=f>tmp?f:tmp;
        b[i][0]=-1;
        i++;
    }
    for(int j=0;j<N;j++){
        A[j+N][j]=-1;
        A[j+2*N][j]=1;
        b[j+2*N][0]=1;
    }
    return f;
}
vector<vector<double> > affine_scaling(vector<vector<double> > &A,vector<vector<double> > &b,vector<vector<double> > &c,int N,int max_iters){
    double gamma=0.5;
    vector<vector<double> > x(N,vector<double>(1,0.5+0.1*rand()/(double)(RAND_MAX)));//0.5+0.1*rand()/(double)(RAND_MAX)
    print_matrix(x,"x.txt");
    for (int i=0;i<max_iters;i++){
        vector<vector<double> > temp=matrix_multi(A,x,3*N,N,1);
        vector<vector<double> > vk=matrix_sub(b,temp,3*N,1);//b-Ax:::::3*N,1
        print_matrix(vk,"vk.txt");
        cout<<"get vk!!!!!"<<endl;
        vector<vector<double> > D=diag(vk,3*N);//3N,3N
        //vector<vector<double> > lambda(3*N,vector<double>(3*N,0));//1e-8
        //for (int i=0;i<3*N;i++) lambda[i][i]=1e-10;
        //D=matrix_plus(D,lambda,3*N,3*N);
        cout<<"get D!!!!!"<<endl;
        vector<vector<double> > INV_D(3*N,vector<double>(3*N));//3N,3N

        if (matrix_inv(D,INV_D,3*N)){
            print_matrix(INV_D,"INV_D.txt");
            vector<vector<double> > INV_D_Square=matrix_multi(INV_D,INV_D,3*N,3*N,3*N);//3N,3N
            vector<vector<double> > trans_A=matrix_trans(A,3*N,N);//N,3N
            vector<vector<double> > matrix_tmp_1=matrix_multi(trans_A,INV_D_Square,N,3*N,3*N);//N,3N
            vector<vector<double> > matrix_tmp_2=matrix_multi(matrix_tmp_1,A,N,3*N,N);//N,N
            vector<vector<double> > INV_matrix_tmp_2(N,vector<double>(N));
            if (matrix_inv(matrix_tmp_2,INV_matrix_tmp_2,N)){
                vector<vector<double> > hx=matrix_multi(INV_matrix_tmp_2,c,N,N,1);//N,1
                print_matrix(hx,"hx.txt");
                vector<vector<double> > hv=matrix_multi(A,hx,3*N,N,1);//3*N,1
                matrix_num_multi(hv,-1,3*N,1);
                print_matrix(hv,"hv.txt");
                if (get_min_value(hv,3*N,1)<0){
                    double alpha=gamma*get_positive_min_value(vk,hv,3*N,1);
                    cout<<"alpha is:"<<alpha<<endl;
                    //double alpha=gamma;
                    matrix_num_multi(hx,alpha,N,1);
                    x=matrix_plus(x,hx,N,3);
                }
                else{
                    cout<<"unbounded!"<<endl;
                    return x;
                }
            }
            else {
                cout<<"failed to inverse"<<endl;
                return x;
            }
        }
        else{
            cout<<"failed to inverse"<<endl;
            return x;
        }
    }
    return x;
}

vector<set<int> > LP_Set_Cover(vector<set<int> > &F,set<int> &X,int N,int max_iters){
    vector<set<int> > res;
    vector<vector<double> > A(3*N,vector<double>(N,0));
    vector<vector<double> > b(3*N,vector<double>(1,0));
    vector<vector<double> > c(N,vector<double>(1,1));
    int f=get_LP_parameters(F,X,A,b,N);
    cout<<endl<<"f:"<<f<<endl;
    print_matrix(A,"A.txt");
    print_matrix(b,"b.txt");
    vector<vector<double> > Xs=affine_scaling(A,b,c,N,max_iters);
    for (int i=0;i<N;i++){
        cout<<Xs[i][0]<<" ";
        if (Xs[i][0]>1.0/(double)f) res.push_back(F[i]);
    }
    return res;
}
*/
set<int> find_min_cost_S(vector<set<int> > &F,set<int> &U_S,vector<double> &price){
    int num=F.size(),index=0;
    set<int> res;
    double minvalue=INT_MAX;
    for (int i=0;i<num;i++){
        int insec=set_sub(F[i],U_S,price);
        if (minvalue>1.0/(double)insec) {
            res=F[i];
            minvalue=1.0/(double)insec;
            index=i;
        }
    }
    F.erase(F.begin()+index);
    addAll(U_S,res);
    return res;
}

vector<set<int> > Dual_LP_Set_Cover(vector<set<int> > F,set<int> X){
    int N=F.size();
    set<int> U_S;
    vector<double> price(N);

    vector<set<int> > res;
    while(U_S.size()!=X.size()){
        set<int> S=find_min_cost_S(F,U_S,price);
        res.push_back(S);
    }
    return res;
}

void Solve_Solution_CH(vector<set<int> > (*func)(vector<set<int> > F,set<int> X),vector<set<int> > &F,set<int> &X,string file){
    clock_t start_time=clock();
    vector<set<int> > C=func(F,X);
    clock_t end_time=clock();
    cout<<"C.size()=="<<C.size()<<endl;
    cout<<"---------------validation-----------"<<endl;
    set<int> U_S;
    for (unsigned int i=0;i<C.size();i++){
        addAll(U_S,C[i]);
    }
    cout<<"U_S.size()=="<<U_S.size()<<endl;
    print_S(U_S,"U_S.txt");
    print_F(C,file.c_str());
    string::size_type idx=file.find('.');
    cout<<"The algorithm of "<<file.substr(0,idx)<<" running time is "<<(double)(end_time-start_time)/(double)(CLOCKS_PER_SEC)<<"s"<<endl;
}

int main()
{
    int N;
    cout<<"Please input the size of X and F:"<<endl;
    cin>>N;
    srand(time(NULL));
    set<int> X=X_generalization(N);//|X|=N,|F|=N
    set<int> Temp=X;
    vector<set<int> > F=F_generalization(N,X);
    vector<set<int> > Temp_F=F;
    X=Temp;
    cout<<"--------------Generalization Sets ends--------------"<<endl;
    print_F(F,"F.txt");
    cout<<"--------------Dual_LP_Set_Cover--------------"<<endl;
    Solve_Solution_CH(Dual_LP_Set_Cover,F,X,"Dual_LP_Set_Cover.txt");
    //cout<<F.size()<<endl;
    F=Temp_F;
    cout<<"--------------Greedy_Set_Cover--------------"<<endl;
    Solve_Solution_CH(Greedy_Set_Cover,F,X,"Greedy_Set_Cover.txt");
    return 0;
}

void print_F(vector<set<int> > F,string file){
    int num=F.size();
    ofstream out(file.c_str());
    for (int i=0;i<num;i++){
        set<int>::iterator S_iter=F[i].begin();
        for (S_iter=F[i].begin();S_iter!=F[i].end();S_iter++){
            out<<*S_iter<<" ";
            //cout<<*S_iter<<" ";
        }
        out<<endl;
        //cout<<endl;
    }
    out.close();
}

void print_S(set<int> S,string file){
    set<int>::iterator S_iter=S.begin();
    ofstream out(file.c_str());
    for (S_iter=S.begin();S_iter!=S.end();S_iter++){
        out<<*S_iter<<" ";
        //cout<<*S_iter<<" ";
    }
    out<<endl;
    out.close();
    //cout<<endl;
}
/*
void print_matrix(vector<vector<double> > &A,string file){
    int m=A.size(),n=A[0].size();
    ofstream out(file.c_str());
    cout<<"------------"<<file<<"-----------:"<<endl;
    for (int i=0;i<m;i++){
        for (int j=0;j<n;j++){
            if(j==0) out<<A[i][j];
            else out<<","<<A[i][j];
            //cout<<A[i][j]<<" ";
        }
        out<<";"<<endl;
        //cout<<endl;
    }
    out.close();
}
*/
