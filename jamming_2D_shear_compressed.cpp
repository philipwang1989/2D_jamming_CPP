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
    // should I use MPI?
    int N;// = 6;
    double alpha;// = 2.0;
    int seed;// = 10;

    N = atoi(argv[1]);
    alpha = atof(argv[2]);
    seed = atoi(argv[3]);

    // step 1. get where compressed states are saved
    string space = "_";
    string savefolder;
    string savename = "jamming_";
    string constPname = "constP_";

    savefolder.append(argv[1]);
    savefolder.append(space);
    savefolder.append(argv[2]);
    savefolder.append(space);
    savefolder.append(argv[3]); // 6_2.0_10
    savename.append(savefolder);// jamming_6_2.0_10
    constPname.append(savefolder); // constP_6_2.0_10
    savename.append(space);     // jamming_6_2.0_10_
    constPname.append(space); // constP_6_2.0_10_
    savefolder.append("/");     // 6_2.0_10/

    int OS = 0;
    string tempname;
    string loadpath;
    string savepath;
    string saveCPpath;
    if (OS == 0) // Linux
    {
        tempname = "/gpfs/loomis/project/ohern/pw374/2D_CPP_FIRE_N_";
        tempname.append(argv[1]);
        tempname.append("_alpha_");
        tempname.append(argv[2]);
        tempname.append("_shearG_scan/");
        loadpath = tempname;
        savepath = tempname;
        loadpath = loadpath.append("compressed_states/"); // ~/compressed_states/
        savepath = savepath.append("sheared_states/"); // ~/sheared_states/
    }
    else // Mac
    {
        tempname = "/Users/philipwang1/Dropbox/Yale/C++/Jamming/";
        loadpath = tempname;
        savepath = tempname;
        loadpath = loadpath.append("compressed_states/"); // ~/compressed_states/
        savepath = savepath.append("sheared_states/"); // ~/sheared_states/
    }

    int p0 = 0;
    bool restart = false; // currently not needed?
    loadpath.append(savefolder); // ~/compressed_states/6_2.0_10/
    savepath.append(savefolder); // ~/sheared_states/6_2.0_10/
    saveCPpath = savepath;
    saveCPpath.append(constPname);
    if (IsPathExist(savepath.c_str()))
    {
        // find all created states and restart from where left behind
        int Nstates;
        if (alpha > 2.0) Nstates = 1600;
        else Nstates = 1000; // alpha = 2.0, 1000 pressure values
        savepath.append(savename); // ~/sheared_states/6_2.0_10/jamming_6_2.0_10_
        for (int p=0; p<Nstates; p++)
        {
            string searchpath = savepath;
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
        mkdir(savepath.c_str(), ACCESSPERMS);
        savepath.append(savename); // ~/sheared_states/6_2.0_10/jamming_6_2.0_10_
    }
    loadpath.append(savename); // ~/compressed_states/6_2.0_10/jamming_6_2.0_10_

    // cout << loadpath << endl;
    // cout << savepath << endl;

    // the real shear part
    // get all Ptarget states, serial    
    int Plow = -7;
    if (alpha > 2.0) Plow = -10;
    int Phigh = -2;

    int Nstates = 200 * (Phigh - Plow);
    double * Plist = new double[Nstates];
    logspace(Plow,Phigh,Nstates,Plist);

    for (int p=p0; p<Nstates; p++)
    {
        double Ptol = Plist[p];
        string loadlocal = loadpath; // ~/compressed_states/6_2.0_10/jamming_6_2.0_10_
        string savelocal = savepath; // ~/sheared_states/6_2.0_10/jamming_6_2.0_10_
        string saveCPlocal = saveCPpath; // ~/sheared_states/6_2.0_10/constP_6_2.0_10_
        stringstream ss;
        ss << p;
        loadlocal.append(ss.str());
        savelocal.append(ss.str());
        saveCPlocal.append(ss.str());
        bool getCPstate = false;
        if (p%20 == 0 || p==Nstates-1) getCPstate = true;
        MD_shearModulus_main(N, alpha, loadlocal.c_str(), savelocal.c_str(), saveCPlocal.c_str(), Ptol, getCPstate);
    }
    return 0;
}