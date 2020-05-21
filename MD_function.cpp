#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <random>
#include <ctime>
#include <sys/stat.h>

using namespace std;

void logspace(int low, int high, int intervals, double * list)
{
    double delta = double(high - low) / double(intervals - 1);
    for (int i=0; i<intervals; i++)
    {
        double curr = low + i*delta;
        list[i] = pow(10, curr);
    }
}

void MD_getRattler(int N, double * Dn, double * x, bool * is_rattler, double gam)
{
    for (int i=0; i<N; i++) is_rattler[i] = false;
    bool if_over = false;
    while (!if_over)
    {
        if_over = true;
        for (int i=0; i<N; i++)
        {
            if (!is_rattler[i])
            {
                vector<double> contact;
                for (int j=0; j<N; j++)
                {
                    if (j!=i && !is_rattler[j])
                    {
                        double sum_r = (Dn[i]+Dn[j])/2.0;
                        double temp_xi = x[2*i];
                        double temp_yi = x[2*i+1];
                        double x_js[4] = {0};
                        double y_js[4] = {0};
                        // virtual 1
                        x_js[0] = x[2*j];
                        y_js[0] = x[2*j+1];
                        // virtual 4
                        if (y_js[0]>temp_yi) y_js[3] = y_js[0] - 1;
                        else y_js[3] = y_js[0] + 1;

                        if (x_js[0]>temp_xi)
                        {
                            if (y_js[0]>temp_yi) x_js[3]=x_js[0]-1-gam;
                            else x_js[3]=x_js[0]-1+gam;
                        }
                        else
                        {
                            if (y_js[0]>temp_yi) x_js[3]=x_js[0]+1-gam;
                            else x_js[3]=x_js[0]+1+gam;
                        }
                        // virtual 2
                        if (y_js[3]<temp_yi) x_js[1] = x_js[3]+gam;
                        else x_js[1] = x_js[3]-gam;
                        y_js[1] = y_js[0];
                        // virtual 3
                        if (y_js[0]>temp_yi) x_js[2] = x_js[0]-gam;
                        else x_js[2] = x_js[0]+gam;
                        y_js[2] = y_js[3];
                        // calculates
                        for (int k=0; k<4; k++)
                        {
                            double temp_xj = x_js[k];
                            double temp_yj = y_js[k];
                            double temp_dy = temp_yj - temp_yi;
                            double temp_dx = temp_xj - temp_xi;
                            if (temp_dx*temp_dx+temp_dy*temp_dy<sum_r*sum_r)
                            {
                                contact.push_back(temp_dx);
                                contact.push_back(temp_dy);
                            }
                        }
                    }
                }
                if (contact.size() > 4)
                {
                    vector<double> theta;
                    for (int k=0; k<contact.size()/2; k++) theta.push_back(atan2(contact[2*k+1],contact[2*k]));
                    sort(theta.begin(),theta.end());
                    if (2*M_PI+theta.front()-theta.back() > M_PI)
                    {
                        is_rattler[i] = true;
                        if_over = false;
                    }
                    for (int k=1; k<theta.size(); k++)
                    {
                        if (theta[k]-theta[k-1]>M_PI)
                        {
                            is_rattler[i] = true;
                            if_over = false;
                        }
                    }
                }
                else
                {
                    is_rattler[i] = true;
                    if_over = false;
                }
                
            }
        }
    }
}

void MD_getStressTensor(const int &N, double * Dn, double * x, double * stress, const double &alpha, const double &gam)
{
    // x runs x0,y0,x1,y1...xN,yN
    int nn, mm;
    double P, dx, dy, Dnm, dnm, F, im, D;
    // double stress[4] = {0}; // Sxx,Sxy,Syx,Syy
    for (int i=0; i<4; i++) stress[i] = 0.0;
    bool * is_rattler = new bool[N];
    for (int i=0; i<N; i++) is_rattler[i] = false;
    MD_getRattler(N, Dn, x, is_rattler, gam);
    D = 1.0;
    for (nn=0; nn<N; nn++)
    {
        D = min(D, Dn[nn]);
        for (mm=nn+1; mm<N; mm++)
        {
            if (!is_rattler[nn] && !is_rattler[mm])
            {
                dy = x[2*mm+1]-x[2*nn+1];
                im = round(dy);
                dy -= im;
                dx = x[2*mm]-x[2*nn];
                dx -= round(dx);
                dx -= (round(dx-im*gam)+im*gam);
                Dnm = 0.5 * (Dn[nn]+Dn[mm]);
                if (fabs(dy)<Dnm)
                {
                    dnm = sqrt(dx*dx+dy*dy);
                    if (dnm < Dnm)
                    {
                        F = -pow((1-dnm/Dnm),alpha-1.0)/Dnm/dnm;
                        stress[0] -= F*dx*dx;
                        stress[1] -= 0.5*F*(dx*dy+dy*dx);
                        stress[2] -= 0.5*F*(dy*dx+dx*dy);
                        stress[3] -= F*dy*dy;
                    }
                }
            }
            
        }
    }
    for (int i=0; i<4; i++) stress[i] *= (D*D/4);
}

double MD_getP(int N, double * Dn, double * x, double alpha, double gam)
{
    // x runs x0,y0,x1,y1...xN,yN
    int nn, mm;
    double P, dx, dy, Dnm, dnm, F, im, D;
    double stress[4] = {0}; // Sxx,Sxy,Syx,Syy
    bool * is_rattler = new bool[N];
    for (int i=0; i<N; i++) is_rattler[i] = false;
    MD_getRattler(N, Dn, x, is_rattler, gam);
    D = 1.0;
    for (nn=0; nn<N; nn++)
    {
        D = min(D, Dn[nn]);
        for (mm=nn+1; mm<N; mm++)
        {
            if (!is_rattler[nn] && !is_rattler[mm])
            {
                dy = x[2*mm+1]-x[2*nn+1];
                im = round(dy);
                dy -= im;
                dx = x[2*mm]-x[2*nn];
                dx -= round(dx);
                dx -= (round(dx-im*gam)+im*gam);
                Dnm = 0.5 * (Dn[nn]+Dn[mm]);
                if (fabs(dy)<Dnm)
                {
                    dnm = sqrt(dx*dx+dy*dy);
                    if (dnm < Dnm)
                    {
                        F = -pow((1-dnm/Dnm),alpha-1.0)/Dnm/dnm;
                        stress[0] -= F*dx*dx;
                        stress[1] -= 0.5*F*(dx*dy+dy*dx);
                        stress[2] -= 0.5*F*(dy*dx+dx*dy);
                        stress[3] -= F*dy*dy;
                    }
                }
            }
            
        }
    }
    for (int i=0; i<4; i++) stress[i] *= (D*D/4);
    P = (stress[0] + stress[3])/2;
    return P;
}

void MD_getCtcNwk(const int &N, double * Dn, double * x, bool * contact_network, const double &gam)
{
    int nn, mm;
    double dx, dy, Dnm, dnm, F, im;
    bool * is_rattler = new bool[N];
    for (int i=0; i<N; i++) is_rattler[i] = false;
    MD_getRattler(N, Dn, x, is_rattler, gam);
    for (int i=0; i<N*N; i++) contact_network[i] = false;
    for (nn=0; nn<N; nn++)
    {
        for (mm=nn+1; mm<N; mm++)
        {
            if (!is_rattler[nn] && !is_rattler[mm])
            {
                dy = x[2*mm+1]-x[2*nn+1];
                im = round(dy);
                dy -= im;
                dx = x[2*mm]-x[2*nn];
                dx -= round(dx);
                dx -= (round(dx-im*gam)+im*gam);
                Dnm = 0.5 * (Dn[nn]+Dn[mm]);
                if (fabs(dy)<Dnm)
                {
                    dnm = sqrt(dx*dx+dy*dy);
                    if (dnm < Dnm) contact_network[nn*N+mm] = true;
                }
            }
        }
    }
}

double MD_getP_DS(const int &N, double * Dn, double * x, const bool * contact_network, const double &alpha, const double &gam)
{
    // x runs x0,y0,x1,y1...xN,yN
    int nn, mm;
    double P, dx, dy, Dnm, dnm, F, im, D;
    double stress[4] = {0}; // Sxx,Sxy,Syx,Syy
    D = 1.0;
    for (nn=0; nn<N; nn++)
    {
        D = min(D, Dn[nn]);
        for (mm=nn+1; mm<N; mm++)
        {
            if (contact_network[nn*N+mm])
            {
                dy = x[2*mm+1]-x[2*nn+1];
                im = round(dy);
                dy -= im;
                dx = x[2*mm]-x[2*nn];
                dx -= round(dx);
                dx -= (round(dx-im*gam)+im*gam);
                Dnm = 0.5 * (Dn[nn]+Dn[mm]);
                dnm = sqrt(dx*dx+dy*dy);
                double delta = 1.0-dnm/Dnm;
                if (delta < 0) // dnm > Dnm - raises ERROR in std::pow
                {
                    F = pow(-delta,alpha-1.0)/Dnm/dnm; // becomes attraction!
                }
                else F = -pow(delta,alpha-1.0)/Dnm/dnm;
                
                if (isnan(F))
                {
                    cout << "NAN F!" << endl;
                    cout << pow((1-dnm/Dnm),alpha-1.0) << endl;
                    cout << (1-dnm/Dnm) << endl;
                    cout << alpha-1.0 << endl;
                } 
                stress[0] -= F*dx*dx;
                stress[1] -= 0.5*F*(dx*dy+dy*dx);
                stress[2] -= 0.5*F*(dy*dx+dx*dy);
                stress[3] -= F*dy*dy;
            }
            
        }
    }
    for (int i=0; i<4; i++) stress[i] *= (D*D/4);
    P = (stress[0] + stress[3])/2;
    return P;
}

double MD_getFtol(const int &N, double * Dn, double * x, double * Fx, const double &alpha, const double &gam)
{
    for (int i=0; i<2*N; i++) Fx[i] = 0.0;
    int nn, mm;
    double Ftol, dx, dy, Dnm, dnm, F, im, D;
    Ftol = 0;
    for (nn=0; nn<N; nn++)
    {
        D = min(D, Dn[nn]);
        for (mm=nn+1; mm<N; mm++)
        {
            dy = x[2 * mm + 1] - x[2 * nn + 1];
            im = round(dy);
            dy -= im;
            dx = x[2 * mm] - x[2 * nn];
            dx -= round(dx);
            dx -= (round(dx - im * gam) + im * gam);
            Dnm = 0.5 * (Dn[nn] + Dn[mm]);
            if (fabs(dy) < Dnm)
            {
                dnm = sqrt(dx * dx + dy * dy);
                if (dnm < Dnm)
                {
                    F = -pow((1 - dnm / Dnm), alpha - 1.0) / Dnm / dnm;
                    Fx[2*nn] += F*dx;
                    Fx[2*mm] -= F*dx;
                    Fx[2*nn+1] += F*dy;
                    Fx[2*mm+1] -= F*dy;
                }
            }
        }
    }    
    for (int i=0; i<2*N; i++) Ftol += (Fx[i] * Fx[i]);
    Ftol /= (N*N);
    return Ftol;
}

double FIRE(int N, double * vx, double * Fx, double a)
{
    // norm - sqrt(sum of the squared)
    double normFx = 0;
    double normFy = 0;
    double normvx = 0;
    double normvy = 0;
    double P = 0;
    for (int i=0; i<2*N; i++) P += vx[i]*Fx[i];

    for (int i=0; i<N; i++)
    {
        normFx += (Fx[2*i]*Fx[2*i]);
        normFy += (Fx[2*i+1]*Fx[2*i+1]);
        normvx += (vx[2*i]*vx[2*i]);
        normvy += (vx[2*i+1]*vx[2*i+1]);
    }
    normFx = sqrt(normFx);
    normFy = sqrt(normFy);
    normvx = sqrt(normvx);
    normvy = sqrt(normvy);

    for (int i=0; i<N; i++)
    {
        vx[2*i] = (1-a) * vx[2*i] + a * (Fx[2*i]/normFx) * normvx;
        vx[2*i+1] = (1-a) * vx[2*i+1] + a * (Fx[2*i+1]/normFy) * normvy;
    }

    return P;
}

void MD_CMzeroing(int N, double * vx, double * m)
{
    double vx_cm, vy_cm, sum_m, Px, Py;
    sum_m = 0.0;
    Px = 0.0;
    Py = 0.0;
    for (int i=0; i<N; i++)
    {
        Px += m[i] * vx[2*i];
        Py += m[i] * vx[2*i+1];
        sum_m += m[i];
    }
    vx_cm = Px / sum_m;
    vy_cm = Py / sum_m;
    for (int i=0; i<N; i++)
    {
        vx[2*i] = vx[2*i] - vx_cm;
        vx[2*i+1] = vx[2*i+1] - vy_cm;
    }
}

double MD_cmpjam_minimization(int N, double * Dn, double * m, double * x, double dt, double alpha, double dphi, double gam, const double& tol)
{
    // first grow particles by dphi
    double totalarea = 0.0;
    for(int i = 0; i < N; i++) totalarea += m[i];
    double rsc = sqrt(1 + dphi/totalarea);
    for (int i = 0; i < N; i++)
    {
        Dn[i] *= rsc;
        m[i] = M_PI * (Dn[i]*Dn[i])/4.0;
    }

    // x runs x0,y0,x1,y1...xN,yN, so as the other arrays with 2N
    long nt, Nt;
    double dt2 = dt * dt;
    
    nt = 0;
    Nt = 5e6 * N/64;
    double * vx = new double[2*N];
    double * ax = new double[2*N];
    double * ax_old = new double[2*N];
    double * Fx = new double[2*N];

    for (int i=0; i<2*N; i++)
    {
        vx[i] = 0.0;
        ax[i] = 0.0;
        ax_old[i] = 0.0;
        Fx[i] = 0.0;
    }

    double P = 0;
    double Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);

    // fire coefficients
    double PP = 0;
    long nfiremin = 5;
    double finc = 1.1;
    double fdec = 0.5;
    double astart = 0.1;
    double a = astart;
    double fa = 0.99;
    double dtmax = 10.0 * dt;
    long cut = nt;

    while (Ftol>tol && nt<Nt)
    {
        // periodic BC, warning!!! does not account for gam != 0.
        
        if (gam == 0)
        {
            for (int i = 0; i < 2*N; i++)
            {
                x[i] = fmod(x[i], 1);
                if (x[i] < 0)
                {
                    x[i] += 1;
                }
            }
        }

        // zeroing center of mass velocity
        if (nt%100 == 0) MD_CMzeroing(N,vx,m);
        
        // first step velocity verlet integration
        for (int i=0; i<2*N; i++) x[i] += vx[i]*dt + ax_old[i]*dt2/2.0;

        // get forces from position
        Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);
        if (Ftol < tol)
        {
            if (nt < 1000) 
            {
                Ftol = tol; // wait until we have a break time
            }
            else
            {
                if (nt >= 1000 && nt < 0.5 * Nt)
                {
                    int Nr = 0;
                    bool * is_rattler = new bool[N];
                    for (int i=0; i<N; i++) is_rattler[i] = false;
                    MD_getRattler(N, Dn, x, is_rattler, gam);
                    Ftol = 0;
                    for (int i=0; i<N; i++) 
                    {
                        if (~is_rattler[i]) 
                        {
                            Nr += 1;
                            Ftol += (Fx[2*i] * Fx[2*i] + Fx[2*i+1] * Fx[2*i+1]);
                        }
                    }
                    Ftol /= (Nr*Nr);
                }
            }
        }

        // calculates acceleration
        for (int i=0; i<N; i++) 
        {
            ax[2*i] = Fx[2*i] / m[i];
            ax[2*i+1] = Fx[2*i+1] / m[i];
        }

        // fire
        if (nt > 0)
        {
            PP = FIRE(N,vx,Fx,a);
            if (PP<0)
            {
                for (int i=0; i<2*N; i++) vx[i] = 0.0;
                cut = nt;
                dt = dt * fdec;
                a = astart;
            }
            else
            {
                if (PP>=0 && nt - cut > nfiremin)
                {
                    dt = min(dt*finc, dtmax);
                    a = a * fa;
                }
            }
            dt2 = dt * dt;
        }

        // second step velocity verlet integration
        for (int i=0; i<2*N; i++) vx[i] += (ax_old[i]+ax[i])*dt/2.0;
        for (int i=0; i<2*N; i++) ax_old[i] = ax[i];
        nt++;
        // if (nt > 1e4 && nt%10000 == 0) cout << nt << ",Ftol=" << Ftol << endl;
    }
    if (nt == Nt) cout << "Nt too small?" << endl;
    P = MD_getP(N,Dn,x,alpha,gam);
    return P;
}

double MD_shrjam_minimization(const int &N, double * Dn, double * m, double * x, double dt, const double &alpha, const double &dG, double &gam, const double& tol)
{
    gam += dG;
    for (int i=0; i<N; i++) x[2*i] += dG * x[2*i+1];
    long nt, Nt;
    double dt2 = dt * dt;
    
    nt = 0;
    Nt = 5e6;
    if (N > 64) Nt = Nt * N/64;

    double * vx = new double[2*N];
    double * ax = new double[2*N];
    double * ax_old = new double[2*N];
    double * Fx = new double[2*N];

    for (int i=0; i<2*N; i++)
    {
        vx[i] = 0.0;
        ax[i] = 0.0;
        ax_old[i] = 0.0;
        Fx[i] = 0.0;
    }
    double P = 0.0;
    double Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);

    if (Ftol<tol) cout << "bad shear!" << endl;

    // fire coefficients
    double PP = 0;
    long nfiremin = 5;
    double finc = 1.1;
    double fdec = 0.5;
    double astart = 0.1;
    double a = astart;
    double fa = 0.99;
    double dtmax = 10.0 * dt;
    long cut = nt;

    while (Ftol>tol && nt<Nt)
    {
        // zeroing center of mass velocity
        MD_CMzeroing(N,vx,m);

        // first step velocity verlet integration
        for (int i=0; i<2*N; i++) x[i] += vx[i]*dt + ax_old[i]*dt2/2.0;

        // get forces from position
        Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);

        // calculates acceleration
        for (int i=0; i<N; i++) 
        {
            ax[2*i] = Fx[2*i] / m[i];
            ax[2*i+1] = Fx[2*i+1] / m[i];
        }

        // fire
        if (nt > 0)
        {
            PP = FIRE(N,vx,Fx,a);
            if (PP<0)
            {
                for (int i=0; i<2*N; i++) vx[i] = 0.0;
                cut = nt;
                dt = dt * fdec;
                a = astart;
            }
            else
            {
                if (PP>=0 && nt - cut > nfiremin)
                {
                    dt = min(dt*finc, dtmax);
                    a = a * fa;
                }
            }
            dt2 = dt * dt;
        }

        // second step velocity verlet integration
        for (int i=0; i<2*N; i++) vx[i] += (ax_old[i]+ax[i])*dt/2.0;
        for (int i=0; i<2*N; i++) ax_old[i] = ax[i];
        nt++;
    }
    if (nt == Nt) cout << "Nt too small?" << endl;
    P = MD_getP(N,Dn,x,alpha,gam);

    delete[] vx;
    delete[] ax;
    delete[] ax_old;
    delete[] Fx;

    return P;
}

void MD_cmpjam_main(int N, double * Dn, double * m, double * x, double dt, double alpha, double Ptol, double dphi, double gam)
{
    double tol = 1e-7;
    double ftol = 1e-30;
    const double tol0 = tol;
    const double ftol0 = ftol;
    const double dphi0 = dphi;
    const double dt0 = dt;

    double P = MD_getP(N,Dn,x,alpha,gam);
    double P0 = P;
    // printf("Current P=%1.5f.\n",P);

    double * Dn_old = new double[N];
    double * m_old = new double[N];
    double * x_old = new double[2*N];
    for (int i=0; i<2*N; i++) x_old[i] = x[i];
    for (int i=0; i<N; i++)
    {
        Dn_old[i] = Dn[i];
        m_old[i] = m[i];
    }

    // save initial state
    double * Dn_ini = new double[N];
    double * m_ini = new double[N];
    double * x_ini = new double[2*N];
    for (int i=0; i<2*N; i++) x_ini[i] = x[i];
    for (int i=0; i<N; i++)
    {
        Dn_ini[i] = Dn[i];
        m_ini[i] = m[i];
    }

    long count = 0;
    long count_max = 5e5;

    bool regular_case = true;
    // bool regular_case = false;

    time_t tstart, tend;
    tstart = time(0);

    while (P < Ptol || P > (1+tol)*Ptol)
    {
        if (P < Ptol)
        {
            dphi = fabs(dphi);
            // copy current to old
            for (int i=0; i<2*N; i++) x_old[i] = x[i];
            for (int i=0; i<N; i++)
            {
                Dn_old[i] = Dn[i];
                m_old[i] = m[i];
            }
            P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
        }
        else if (P > (1+tol)*Ptol)
        {
            if (regular_case) // regular case
            {
                // copy old to current
                for (int i=0; i<2*N; i++) x[i] = x_old[i];
                for (int i=0; i<N; i++)
                {
                    Dn[i] = Dn_old[i];
                    m[i] = m_old[i];
                }
                dphi /= 2.0;
            }
            else // difficult case
            {
                mt19937 generator(time(0));
                uniform_real_distribution<double> distribution(0.0, 1.0);
                double randscale;
                do {randscale = distribution(generator);} while (randscale < 0.1);
                dphi = -fabs(dphi) * randscale;
            }            
            
            P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
        }
        // print out here
        double totalarea = 0.0;
        for(int i = 0; i < N; i++) totalarea += m[i];

        /*
        if (totalarea > 0.8 && count > 0 && count < 150) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        else if (totalarea > 0.8 && count > 150 && count < 1e3 && count%100 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        
        else if (totalarea > 0.8 && count > 1e3 && count < 1e4 && count%1000 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        else if (totalarea > 0.8 && count > 1e4 && count%10000 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        */

        count ++;
      
        if (count > count_max || fabs(dphi) < 3.0e-16)
        {   
            count = 0;
            regular_case = false;

            mt19937 generator(time(0));
            uniform_real_distribution<double> distribution(0.0, 1.0);
            double randscale;
            do {randscale = distribution(generator);} while (randscale < 0.1);
            dphi = -dphi0 * randscale;
            P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);

            double totalarea = 0.0;
            for(int i = 0; i < N; i++) totalarea += m[i];
            // printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        }
    }
    double totalarea = 0.0;
    for(int i = 0; i < N; i++) totalarea += m[i];
    tend = time(0);
    printf("%ld,t=%1.1fs,phi=%1.7f,dphi=%e,P=%e.\n",count,difftime(tend,tstart),totalarea,dphi,P);
}

void MD_cmpjam_main_DS(int N, double * Dn, double * m, double * x, bool * contact_network, double dt, double alpha, double Ptol, double dphi, double gam)
{
    double tol = 1e-7;
    double ftol = 1e-30;
    const double tol0 = tol;
    const double ftol0 = ftol;
    const double dphi0 = dphi;
    const double dt0 = dt;

    double P = MD_getP_DS(N,Dn,x,contact_network,alpha,gam);
    double P0 = P;
    // printf("Initial P=%e.\n",P);

    double * Dn_old = new double[N];
    double * m_old = new double[N];
    double * x_old = new double[2*N];
    for (int i=0; i<2*N; i++) x_old[i] = x[i];
    for (int i=0; i<N; i++)
    {
        Dn_old[i] = Dn[i];
        m_old[i] = m[i];
    }

    // save initial state
    double * Dn_ini = new double[N];
    double * m_ini = new double[N];
    double * x_ini = new double[2*N];
    for (int i=0; i<2*N; i++) x_ini[i] = x[i];
    for (int i=0; i<N; i++)
    {
        Dn_ini[i] = Dn[i];
        m_ini[i] = m[i];
    }

    long count = 0;
    long count_max = 5e5;

    bool regular_case;

    if (P > (1+tol)*Ptol) regular_case = false; // loaded a overcompressed state
    else regular_case = true;

    if (!regular_case) cout << "Overcompressed initial state!" << endl;

    time_t tstart, tend;
    tstart = time(0);

    while (P < Ptol || P > (1+tol)*Ptol)
    {
        if (P < Ptol)
        {
            dphi = fabs(dphi);
            // copy current to old
            for (int i=0; i<2*N; i++) x_old[i] = x[i];
            for (int i=0; i<N; i++)
            {
                Dn_old[i] = Dn[i];
                m_old[i] = m[i];
            }
            P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
        }
        else if (P > (1+tol)*Ptol)
        {
            if (regular_case) // regular case
            {
                // copy old to current
                for (int i=0; i<2*N; i++) x[i] = x_old[i];
                for (int i=0; i<N; i++)
                {
                    Dn[i] = Dn_old[i];
                    m[i] = m_old[i];
                }
                dphi /= 2.0;
            }
            else // difficult case
            {
                mt19937 generator(time(0));
                uniform_real_distribution<double> distribution(0.0, 1.0);
                double randscale;
                do {randscale = distribution(generator);} while (randscale < 0.1);
                dphi = -fabs(dphi) * randscale;
            }            
            
            P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);
        }
        // update P here based on contact network
        P = MD_getP_DS(N,Dn,x,contact_network,alpha,gam);

        // print out here
        double totalarea = 0.0;
        for(int i = 0; i < N; i++) totalarea += m[i];

        /*
        if (totalarea > 0.8 && count > 50 && count < 150) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        else if (totalarea > 0.8 && count > 150 && count < 1e3 && count%100 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        else if (totalarea > 0.8 && count > 1e3 && count < 1e4 && count%1000 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        else if (totalarea > 0.8 && count > 1e4 && count%10000 == 0) printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        */

        count ++;
      
        if (count > count_max || fabs(dphi) < 3.0e-16)
        {   
            count = 0;
            regular_case = false;

            mt19937 generator(time(0));
            uniform_real_distribution<double> distribution(0.0, 1.0);
            double randscale;
            do {randscale = distribution(generator);} while (randscale < 0.1);
            dphi = -dphi0 * randscale;
            P = MD_cmpjam_minimization(N,Dn,m,x,dt,alpha,dphi,gam,ftol);

            double totalarea = 0.0;
            for(int i = 0; i < N; i++) totalarea += m[i];
            // printf("%ld,phi=%1.7f,dphi=%e,P=%e.\n",count,totalarea,dphi,P);
        }
    }
    double totalarea = 0.0;
    for(int i = 0; i < N; i++) totalarea += m[i];
    tend = time(0);
    printf("%ld,t=%1.1fs,phi=%1.7f,dphi=%e,P=%e.\n",count,difftime(tend,tstart),totalarea,dphi,P);
}

void scale(int N, double * Dn, double * m, double rescale)
{
    double D2;
    for(int i = 0; i < N; i++)
    {
        Dn[i] = Dn[i] * rescale;
        D2 = Dn[i] * Dn[i];
        m[i] = M_PI * D2/ 4.0;
    }
}

//check if any particles are touching
bool anytouch(int N, double * pos, double * sc)
{
	for(int i=0; i<N-1; i++)
	{
		for(int j=i+1; j<N; j++)
		{
			//vector pointing from particle i to j
			double rij[2];
			//periodic boundary conditions
			rij[0] = pos[2*j]-pos[2*i]-round(pos[2*j]-pos[2*i]);
			rij[1] = pos[2*j+1]-pos[2*i+1]-round(pos[2*j+1]-pos[2*i+1]);
			//sum of 2 radii
			double bigR = sc[i]+sc[j];
			//they touch if the distance between them is less than bigR
			if(rij[0]*rij[0]+rij[1]*rij[1] < bigR*bigR) return true;
		}
	}
	return false;
}

// check if file exists
bool fileExist(const char * name)
{
    if (FILE *file = fopen(name, "r")) 
    {
        fclose(file);
        return true;
    } 
    else 
    {
        return false;
    }
}

// check if path exists
bool IsPathExist(const char * s)
{
  struct stat buffer;
  return (stat (s, &buffer) == 0);
}

// write to output
void writeResult(const char * path, int N, double * Dn, double * m, double * x)
{
    FILE *xyz = fopen(path, "w+");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", Dn[j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", m[j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", x[2*j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", x[2*j+1]);
    fprintf(xyz, "\n");

    fclose(xyz);
}

// write to output
void writeResultToFile(FILE * xyz, const int &N, double * Dn, double * m, double * x)
{
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", Dn[j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", m[j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", x[2*j]);
    fprintf(xyz, "\n");
    for (int j=0; j<N; j++) fprintf(xyz, "%1.20f ", x[2*j+1]);
    fprintf(xyz, "\n");
}

// load result
void loadResult(const char * path, int N, double * Dn, double * m, double * x)
{
    FILE *xyz = fopen(path, "r");
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &Dn[j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &m[j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[2*j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[2*j+1]);
    fclose(xyz);
}

void loadResultFromFile(FILE * xyz, const int &N, double * Dn, double * m, double * x)
{
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &Dn[j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &m[j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[2*j]);
    for (int j=0; j<N; j++) fscanf(xyz, "%lf", &x[2*j+1]);
}

void MD_shearModulus_main(const int &N, const double &alpha, const char * loadpath, const char * savepath, const char * saveCPpath, const double &Ptol, const bool &getCPstate, const bool &positive_shear)
{
    FILE *abc = fopen(loadpath, "r");
    if (abc == NULL)
    {
        cout << loadpath << " not exists!" << endl;
        return;
    }
    // array
    double * Dn = new double[N];
    double * m = new double[N];
    double * x = new double[2*N];
    loadResultFromFile(abc, N, Dn, m, x);
    fclose(abc);

    // measurements
    double * Fx = new double[2*N];
    double stress[4] = {0};

    // simulation parameters
    double dG;
    int N_shear_steps = 20;
    int curr_steps = 0;
    double dt = 5e-4;
    double gam = 0.0;
    double ftol = 1e-30;

    if (positive_shear) dG = 1e-9;
    else dG = -1e-9;

    // open file to save sheared state, write gam=0 state
    FILE *xyz = fopen(savepath, "w+");
    writeResultToFile(xyz, N, Dn, m, x);
    FILE *cpstate = NULL;
    if (getCPstate)
    {
        cpstate = fopen(saveCPpath, "w+");
        writeResultToFile(cpstate, N, Dn, m, x);
    }

    while (curr_steps < N_shear_steps)
    {
        // shear at constant volume to gam+dG
        MD_shrjam_minimization(N, Dn, m, x, dt, alpha, dG, gam, ftol);
        for (int i=0; i<2*N; i++) Fx[i] = 0.0;
        double Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);
        MD_getStressTensor(N, Dn, x, stress, alpha, gam);
        cout << "Ftol=" << Ftol << ",sigma_xy=" << -stress[1] << endl;

        // write to output file, do not close
        writeResultToFile(xyz, N, Dn, m, x);

        // find fixed pressure state if necessary
        // before compression need to save the state
        if (getCPstate)
        {
            // save initial state
            double * Dn_ini = new double[N];
            double * m_ini = new double[N];
            double * x_ini = new double[2*N];
            for (int i=0; i<2*N; i++) x_ini[i] = x[i];
            for (int i=0; i<N; i++)
            {
                Dn_ini[i] = Dn[i];
                m_ini[i] = m[i];
            }
            double dphi = 1e-5;
            MD_cmpjam_main(N, Dn, m, x, dt, alpha, Ptol, dphi, gam);

            // write fixed pressure state
            writeResultToFile(cpstate, N, Dn, m, x);
            
            // recover initial sheared state
            for (int i=0; i<2*N; i++) x[i] = x_ini[i];
            for (int i=0; i<N; i++)
            {
                Dn[i] = Dn_ini[i];
                m[i] = m_ini[i];
            }

            delete[] Dn_ini;
            delete[] m_ini;
            delete[] x_ini;
        }      
        cout << "Step=" << curr_steps << ",gamma=" << gam << ",P=" << Ptol << " completed." << endl;
        curr_steps += 1;
    }

    fclose(xyz);
    if (cpstate != NULL) fclose(cpstate);

    delete[] Dn;
    delete[] m;
    delete[] x;
    delete[] Fx;
}

void MD_mapping_shear_func(const int &N, const double &alpha, const char * loadpath, double * Plist, const double &dG)
{
    FILE *abc = fopen(loadpath, "r");
    if (abc == NULL)
    {
        cout << loadpath << " not exists!" << endl;
        return;
    }
    // array
    double * Dn = new double[N];
    double * m = new double[N];
    double * x = new double[2*N];
    loadResultFromFile(abc, N, Dn, m, x);
    fclose(abc);

    // measurements
    double * Fx = new double[2*N];
    double stress[4] = {0};

    // step 1. shear at const. volume, to get P_list
    // always load p=0 state -> correspond to lowest pressure

    // simulation parameters
    // double dG = 1e-9;
    int N_shear_steps = 20;
    int curr_steps = 0;
    double dt = 5e-4;
    double gam = 0.0;
    double ftol = 1e-30;
    int Nstates = 21;

    // get contact network
    bool * contact_network = new bool[N*N];
    MD_getCtcNwk(N, Dn, x, contact_network, gam);

    Plist[0] = MD_getP(N, Dn, x, alpha, gam);

    while (curr_steps < N_shear_steps)
    {
        // shear at constant volume to gam+dG
        MD_shrjam_minimization(N, Dn, m, x, dt, alpha, dG, gam, ftol);
        for (int i=0; i<2*N; i++) Fx[i] = 0.0;
        double Ftol = MD_getFtol(N,Dn,x,Fx,alpha,gam);
        MD_getStressTensor(N, Dn, x, stress, alpha, gam);
        cout << "Ftol=" << Ftol << ",sigma_xy=" << -stress[1] << endl;
        
        curr_steps += 1;
        Plist[curr_steps] = MD_getP_DS(N, Dn, x, contact_network, alpha, gam);
        // cout << "Step=" << curr_steps << ",gamma=" << gam << ",Pds=" << Plist[curr_steps] << ",P=" << MD_getP(N, Dn, x, alpha, gam) << " completed." << endl;
        // cout << "Step=" << curr_steps << ",gamma=" << gam << ",P=" << Plist[curr_steps] << " completed." << endl;
    }

    delete[] Dn;
    delete[] m;
    delete[] x;
    delete[] Fx;

    cout << "Done shearing!" << endl;
}

void MD_mapping_shearCP_func(const int &N, const double &alpha, const char * loadpath, const char * savepath, const double &Ptol, const double &dG)
{
    FILE *abc = fopen(loadpath, "r");
    if (abc == NULL)
    {
        cout << loadpath << " not exists!" << endl;
        return;
    }
    // array
    double * Dn = new double[N];
    double * m = new double[N];
    double * x = new double[2*N];
    loadResultFromFile(abc, N, Dn, m, x);
    fclose(abc);

    // Simulation parameters
    double dt = 5e-4;
    double gam = 0.0;
    double ftol = 1e-30;

    // get contact network
    bool * contact_network = new bool[N*N];
    MD_getCtcNwk(N, Dn, x, contact_network, gam);

    // step 1: compress to target pressure and write to file, at zero gamma
    double dphi = 1e-5;
    MD_cmpjam_main_DS(N, Dn, m, x, contact_network, dt, alpha, Ptol, dphi, gam);
    FILE *xyz = fopen(savepath, "w+");
    writeResultToFile(xyz, N, Dn, m, x);

    // step 2: shear and find const. pressure state
    // double dG = 1e-9;
    int N_shear_steps = 20;
    int curr_steps = 0;
    while (curr_steps < N_shear_steps)
    {
        // shear at constant volume to gam+dG
        MD_shrjam_minimization(N, Dn, m, x, dt, alpha, dG, gam, ftol);

        // get const. P state
        // MD_cmpjam_main(N, Dn, m, x, dt, alpha, Ptol, dphi, gam);
        MD_cmpjam_main_DS(N, Dn, m, x, contact_network, dt, alpha, Ptol, dphi, gam);

        // write to file
        writeResultToFile(xyz, N, Dn, m, x);

        cout << "Step=" << curr_steps << ",gamma=" << gam << ",P=" << MD_getP(N, Dn, x, alpha, gam) << " completed." << endl;
        curr_steps += 1;
    }

    fclose(xyz);

    delete[] contact_network;
}