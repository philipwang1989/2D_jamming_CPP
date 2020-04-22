#include <iostream>
#include <stdio.h>
#include <cstring>
#include <random>
#include <math.h>
#include <algorithm>
#include <sys/stat.h>
#include <sstream>
#include "MD_function.hpp"

using namespace std;

int main(int argc, char ** argv)
{
    /*
    to compile on cluster:
    g++ jamming_2D_prep_comp.cpp MD_function.hpp MD_function.cpp -static -O3 -o jamming.out
    local machine:
    clang++ jamming_2D_prep_comp.cpp MD_function.hpp MD_function.cpp -O3
    */
    int N;// = 6;
    double alpha;// = 1.0;
    int seed;// = 10;
    
    N = atoi(argv[1]);
    alpha = atof(argv[2]);
    seed = atoi(argv[3]);

    //Save output file
    string space = "_";
    string savepath;
    string savefolder;
    string savename = "jamming_";
    string restartfile; // for alpha>2 restarting
    // string filenamexyz;

    savefolder.append(argv[1]);
    savefolder.append(space);
    savefolder.append(argv[2]);
    savefolder.append(space);
    savefolder.append(argv[3]);
    savename.append(savefolder);
    savename.append(space);
    savefolder.append("/");

    int OS = 0;
    string filename;
    if (OS == 0) // Linux
    {
        filename = "/gpfs/loomis/project/ohern/pw374/2D_CPP_FIRE_N_";
        filename.append(argv[1]);
        filename.append("_alpha_");
        filename.append(argv[2]);
        filename.append("_shearG_scan/compressed_states/");
    }
    else // Mac
    {
        filename = "/Users/philipwang1/Dropbox/Yale/C++/Jamming/";
    }

    // get restart file name and location
    if (alpha > 2.0)
    {
        if (OS == 0)
        {
            restartfile = "/gpfs/loomis/project/ohern/pw374/2D_CPP_FIRE_N_";
            restartfile.append(argv[1]);
            restartfile.append("_alpha_2.0");
            restartfile.append("_shearG_scan/compressed_states/");
        }
        else
        {
            restartfile = filename;
        }
        string loadfolder = savefolder;
        size_t found = loadfolder.find("2.5");
        loadfolder.replace(found, 3, "2.0");
        restartfile.append(loadfolder);
        string loadname = savename;
        found = loadname.find("2.5");
        loadname.replace(found, 3, "2.0");
        loadname.append("0");
        restartfile.append(loadname);
        // cout << restartfile << endl;
    }

    int p0 = 0;
    bool restart = false;
    filename.append(savefolder);
    if (IsPathExist(filename.c_str()))
    {
        // find all created states and restart from where left behind
        int Nstates;
        if (alpha > 2.0) Nstates = 1600;
        else Nstates = 1000; // alpha = 2.0, 1000 pressure values
        filename.append(savename);

        for (int p=0; p<Nstates; p++)
        {
            string searchpath = filename;
            stringstream ss;
            ss << p;
            searchpath.append(ss.str());
            if (p == Nstates-1 && fileExist(searchpath.c_str()))
            {
                cout << searchpath << " completed!" << endl;
                return 0;
            }
            if (!fileExist(searchpath.c_str())) // file does not exist
            {
                p0 = p-2;
                restart = true;
                if (p0 < 0) 
                {
                    p0 = 0;
                    restart = false;
                }
                cout << "Restart with " << searchpath << endl;
                break;
            }
        }
    }
    else 
    {
        mkdir(filename.c_str(), ACCESSPERMS);
        filename.append(savename);
    }
    // filenamexyz = filename;

	char *filechar = new char[filename.length() + 1];
	strcpy(filechar, filename.c_str());
    printf("Output char is: %s\n", filechar);

	// char *filexyzchar = new char[filenamexyz.length() + 1];
	// strcpy(filexyzchar, filenamexyz.c_str());	
    // printf("XYZ char is: %s\n", filexyzchar);

    mt19937 generator(seed);
    uniform_real_distribution<double> distribution(0.0, 1.0); // Generate number between 0, 1

    // array
    double * Dn = new double[N];
    double * m = new double[N];
    double * x = new double[2*N];

    double D = 1.0;
    double G = 1.4;
    double gam = 0.0;

    double Lx = 1.0;
    double Ly = 1.0;

    // initial values
    for (int i=0; i<N/2; i++) Dn[i] = D;
    for (int i=N/2; i<N; i++) Dn[i] = G;

    // rescale particles to ipf
    double ipf = 0.01;
    double totalarea = 0.0;
    for(int i = 0; i < N; i++) totalarea += M_PI * Dn[i] * Dn[i] / 4.0;
    double phitot = totalarea / (Lx * Ly);
    double rescale = sqrt(ipf / totalarea);
    scale(N, Dn, m, rescale);
    totalarea = 0.0;
    for(int i = 0; i < N; i++) totalarea += m[i];
    phitot = totalarea / (Lx * Ly);
    // printf("Iinitial phi = %5.10f\n", phitot);

    // simulation parameters
    double dt = 5e-4;
    double Ptol = 1e-7;
    double dphi = 1e-2;

    if (alpha > 2.0) Ptol = 1e-10;

    // get all Ptarget states, serial    
    int Plow = -7;
    if (alpha > 2.0) Plow = -10;
    int Phigh = -2;

    int Nstates = 200 * (Phigh - Plow);
    double * Plist = new double[Nstates];
    logspace(Plow,Phigh,Nstates,Plist);

    if (!restart)
    {
        // Initial particle position
        cout << "Fresh new start!" << endl;
        cout << filename << endl;
        if (alpha < 2.5) // for alpha=2.0
        {
            do
            {
                for (int i=0; i<2*N; i++) x[i] = distribution(generator);
            } while (anytouch(N, x, Dn));
            // x runs x0,y0,x1,y1...xN,yN
        }
        else
        {
            // step 1: check if alpha=2.0, zero pressure state exists
            // step 2: if exist, load the state, else return error message
            if (fileExist(restartfile.c_str()))
            {
                loadResult(restartfile.c_str(),N,Dn,m,x);
                for(int i = 0; i < N; i++) totalarea += m[i];
                cout << "Loaded from alpha=2.0! Current P=" << MD_getP(N, Dn, x, alpha, 0.0) << ",phi=" << totalarea << "!" << endl;
                if (MD_getP(N, Dn, x, alpha, 0.0) < Ptol) dphi = 1e-10;
                else dphi = -1e-4;
            }
            else
            {
                cout << "No available 0 pressure state!" << endl;
                return 0;
            }            
        }        
        // compress to approx. 0 pressure
        MD_cmpjam_main(N, Dn, m, x, dt, alpha, Ptol, dphi, gam);
    }
    else
    {
        // reload the previous acceptable configuration
        cout << "Restarting and reloading!" << endl;
        string load = filename;
        stringstream ss;
        if (p0-1 >= 0) ss << (p0-1);
        else ss << p0;        
        load.append(ss.str());
        loadResult(load.c_str(),N,Dn,m,x);
        cout << "Reloaded!" << endl;
    }

    /*
    for (int j=0; j<N; j++) printf("%1.16f ", Dn[j]);
    printf("\n");
    for (int j=0; j<N; j++) printf("%1.16f ", x[2*j]);
    printf("\n");
    for (int j=0; j<N; j++) printf("%1.16f ", x[2*j+1]);
    printf("\n");
    */

    dphi = 1e-5;
    for (int p=p0; p<Nstates; p++)
    {
        Ptol = Plist[p];
        MD_cmpjam_main(N, Dn, m, x, dt, alpha, Ptol, dphi, gam);
        // write to file here
        string local = filename;
        stringstream ss;
        ss << p;
        local.append(ss.str());
        writeResult(local.c_str(),N,Dn,m,x);
    }
    return 0;
}