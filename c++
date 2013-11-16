#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
using namespace arma;
using namespace std;

double pi = atan(1)*4;

double least_squares(mat exact, mat numerical,int n){
    double difference=0;
    double d;
    for (int i=0;i<n;i++){
        d = exact(0,i)-numerical(0,i);
        difference += d*d;
    }
    return difference;
}

mat analytical(mat x, double t,int m,double k){
    double An;
    double kk=k*k;
    mat v =zeros<mat>(1,m);

    for (int i=30; i>0; i--){
        double n= i;
        An = -2.0/n/pi;
        v += An*sin(pi*n*x)*exp(-n*n*kk*t);
    }
    return v;
}

mat tridiag(double sub, double diag, double super, int n, mat v){

    vec a=vec(n);
    vec b=vec(n);
    vec c=vec(n);

    //setter verdier inn i vektorene over:
    for (int i=0;i<n;i++){
        a(i)=sub;
        b(i)=diag;
        c(i)=super;
    }

    double btemp = b(0);
    vec temp = vec(n);
    mat u=mat(1,n);
    u(0,0) = v(0,0)/btemp;

    for(int i=1 ; i < n ; i++) {
        temp(i) = c(i-1)/btemp;
        btemp = b(i)-a(i)*temp(i);
        u(0,i) = (v(i) - a(i)*u(0,i-1))/btemp;
    }

    // bakoversubstitusjon
    for(int i=n-2 ; i >= 0 ; i--) {
        u(0,i) -= temp(i+1)*u(0,i+1);
        }

        return u;
    }

mat f_forward_euler(double dx, double dt, int m, mat v,double rand0, double rand_end){
    double alpha=dt/dx/dx;
    mat vnew= zeros<mat>(1,m);
    for (int i=1;i<m-1;i++){
        vnew(0,i) = alpha * v(0,i-1) + (1 - 2*alpha) * v(0,i) + alpha * v(0,i+1);
    }
    vnew(0,0)=rand0;
    vnew(0,m-1)=rand_end;
    return vnew;
}

mat f_backward_euler(double dx, double dt, int n, mat v,
                     double rand0, double rand_end){

    double alpha=dt/dx/dx;
    double diagonal=1+2*alpha;
    double subdiagonal=-alpha;
    double superdiagonal=-alpha;

    v = tridiag(subdiagonal,diagonal,superdiagonal,n,v);
    v(0,0)=rand0;
    v(0,n-1)=rand_end;
    return v;
}

mat f_crank_nicolson(double dx, double dt, int n, mat u, double rand0, double rand_end){
    //koden er hentet fra laereboka, men tilpasset programmet mitt
    double alpha=dt/dx/dx;
    double a=-alpha;
    double c=-alpha;
    double b = 2 + 2*alpha;

    mat r=mat(1,n);

    for (int i = 1; i < n-1; i++) {
        r(0,i) = alpha*u(0,i-1) + (2 - 2*alpha)*u(0,i) + alpha*u(0,i+1);
    }

    r(0,0) = 0;
    r(0,n-1) = 0;
    // Then solve the tridiagonal matrix
    u=tridiag(a,b,c,n,r);
    u(0,0) = rand0;
    u(0,n-1) = rand_end;
    return u;
}

///////////
//   MAIN:
////////////
int main()
{
    cout.precision(16);
    // --- initialisering:    

    double dx=0.01;
    double dt=5.0e-5;
    //double dx=0.001;
    //double dt=0.001;

    double T=1.0;

    double L=1;
    double k=pi;
    double rand_max=0.0;
    double rand_min=0.0;

    int n= T/dt+1;
    int m= L/dx+1;

    mat x= mat(1,m);
    for(int i=0;i<m;i++){
        x(0,i)=i*dx;
    }
    mat t= mat(1,n);
    for(int i=0;i<n;i++){
        t(0,i)=i*dt;
    }

    // ***** initialiserer v(x,0) ****
    mat v=zeros<mat>(1,m);
    v(0,0)=1.0; // u(0,0)=1.0
    v += -1+x;    // v(x,0)= u(x,0)-(1-x)
    // *******************************

    mat v_forward_euler  = v;
    mat v_backward_euler = v;
    mat v_crank_nicolson = v;

    string filnavn="analytisk.txt";
    string filnavn2="forward.txt";
    string filnavn3="backward.txt";
    string filnavn4="crank_nicolson.txt";
    string filnavn5="timearray.txt";

    ofstream analytisk, forward, backward, crank_nicolson,timearray;

    analytisk.open (filnavn.c_str());
    forward.open (filnavn2.c_str());
    backward.open (filnavn3.c_str());
    crank_nicolson.open (filnavn4.c_str());
    timearray.open (filnavn5.c_str());

    analytisk      << n<<"  "<<m << "  "<<T<<"  "<<dt<<"  "<<L<<"  "<<dx <<endl;
    forward        << n<<"  "<<m << "  "<<T<<"  "<<dt<<"  "<<L<<"  "<<dx <<endl;
    backward       << n<<"  "<<m << "  "<<T<<"  "<<dt<<"  "<<L<<"  "<<dx <<endl;
    crank_nicolson << n<<"  "<<m << "  "<<T<<"  "<<dt<<"  "<<L<<"  "<<dx <<endl;
    // ------ iterering: -----


    analytisk << analytical(x,t(0),m,k);
    forward   << v;
    backward  << v;
    crank_nicolson << v;
    timearray << t(0,0)<<endl;
    int counter=1;

    for(int i=1;i<n;i++){ // n tidsteg
        v_forward_euler  = f_forward_euler(dx,dt,m,v_forward_euler,rand_min,rand_max);
        v_backward_euler = f_backward_euler(dx,dt,m,v_backward_euler,rand_min,rand_max);
        v_crank_nicolson = f_crank_nicolson(dx,dt,m,v_crank_nicolson,rand_min,rand_max);

        if(counter%1000==0){
        analytisk << setw(17) << analytical(x,t(i),m,k);
        forward        << setw(17) << v_forward_euler;

        backward       << setw(17) << v_backward_euler;
        crank_nicolson << setw(17) << v_crank_nicolson;
        timearray<< t(0,i)<<endl;
        }
        counter++;
    }
    // ----------------------
    return 0;
}
