#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <math.h>
#include <iostream>
#include "mpi.h"
#include "myfftw.h"
#include "file_handling.h"
#include "lmp_handling.h"

class crystal{
    int id;
    char el[4];
    double coord[3];
};

class mol{
    int natom;
    double cube[6];
    double *coord;
};

void extract_velocities_each_atom(
        char* filename, double* velocities,
        int nmd, int natoms,
        int atom_index = 0, char* label2read = "vx", const int nskip=9);

void cal_Irange_MPI(int *Iinit, int *Ifin, int *neach, int myrank, int nprocs, int Nloop);

void mkinput_sed(char *file);

//--- for read dump file
//void get_ROT_INFO(int ia, int nuat, int Nrot, int *irot, int *isite);
//void get_ROT_INFO2(int iatom, int nuat, int Nrot, int *iz, int *irot, int *isite);
int get_irow2read(const std::string &str1, const std::string &sword);
double extract_data_from_line(const std::string &str1, int Iget);
//void extract_isite_velo_for_fft3D(char *FDUMP, int ISITE_GET, int SEDTYPE, int ITYPE_CNT, 
//        int Nall, int Nz, int Nrot, int Nmd, int Nsite, double *velo_site, double radius);
void conv_xyz2rtz(int nat, double **xyz0, double **rtz0, double *center);
void substract_cylind(int Nat, double **vrtz, double **rtz0, double **rtz1);

//--- results
void output_sed(char *krotfile, char *sumfile, char *LABEL, int nk, int nc, int nw, double *wlist, double *sed);

int main(int argc, char **argv){
    
    //---- start MPI
    int nprocs, myrank;
    int Iinit, Ifin, neach;
    //MPI_Init(&argc, &argv);
    //MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    //MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    //if(myrank==0){
    printf("---------------------------------\n");
    printf("\n");
    printf("       dump2sed ver.1.0\n");
    printf("\n");
    printf("---------------------------------\n");
    //}
    
    //--- parameters
    int i, j, ii, ia, idat, count;
    int Nuat, Nrot, i_coord_type;
    double tstep;
    char IFILE[128], word[128];
    char FDUMP[128], FXYZ[128], CTYPE[128];
    char line[256];

    //--- input file name
    if (argc > 1) sprintf(IFILE, "%s", argv[1]);
    else{
        sprintf(IFILE, "example_sed.in");
        printf("Input file name.\nSee %s\n", IFILE);
        mkinput_sed(IFILE);
        exit(0);
    }
    
    //--- file check
    ifstream ifs_check(IFILE);
    if(ifs_check.fail()){
        printf("Cannot find %s\n", IFILE);
        exit(0);
    }
    ifs_check.close();
    
    //--- read input file
    int ised;
    file::searchword(IFILE, "DUMP_FILE", FDUMP);
    file::searchword(IFILE, "INIT_COORD", FXYZ);
    file::searchword(IFILE, "TIMESTEP", line);   tstep = atof(line);
    int ITYPE_CNT = 1;
    
    //--- read initial xyz file
    int Natom = lmp::get_natom_from_dump(FDUMP);
    int Ntype = lmp::get_ntype_from_input(FXYZ);
    int *ID_all = new int[Natom];        //--------- new ID_ALL[Natom]
    int *TY_all = new int[Natom];        //--------- new TY_all[Natom]
    double cube[6];
    double *mass = new double[Ntype];    //--------- new mass[Ntype]
    double **Call = new double *[Natom];
    for(i=0;i<Natom;i++) Call[i] = new double[3];
    lmp::get_coord_from_input(FXYZ, Ntype, mass, Natom, ID_all, TY_all, cube, Call);
    
    //--- read dump file
    int nstep, Nmd, nline = 0;
    nstep = lmp::get_dumpstep(FDUMP);
    nline = file::get_num_of_lines(FDUMP);
    Nmd   = nline / (Natom + 9);
    Nmd   = pow(2, int(log2(Nmd)));    //--- Nmd = 2^**
    
    // --- prepare arrays
    double *velo_each = new double[Nmd];
    double *dos_each = new double[Nmd];
    double *dos_total = new double[Nmd];
    char labels[3][10] = {"vx", "vy", "vz"};
    
    for(i=0; i<Nmd; i++)
        dos_total[i] = 0.;
    
    for(ia=0; ia<Natom; ia++){
        for(j=0; j<3; j++){
            
            cout << ia << " " << j << endl;
            
            extract_velocities_each_atom(FDUMP, velo_each, Nmd, Natom, ia, labels[j]);
            
            MYFFTW::fftw1d(Nmd, velo_each, dos_each);
            
            for(idat=0; idat<Nmd; idat++)
                dos_total[idat] += velo_each[idat] / double(Natom * 3.);
            
            //for (int imd=0; imd<Nmd; imd++)
            //    //cout << velo_each[imd] << " " << dos_each[imd] << endl;
            //    cout << dos_each[imd] << endl;
        }
    }
     
    //delete[] velo;
    MPI_Finalize();
    return 0;
}

void cal_Irange_MPI(int *Iinit, int *Ifin, int *neach, int myrank, int nprocs, int Nloop){
    int iproc;
    int I0, I1, ncal;
    I0 = I1 = 0;
    for(iproc=0; iproc<myrank+1; iproc++){
        ncal = int(Nloop / nprocs);
        if(iproc < Nloop - int(Nloop/nprocs)*nprocs) ncal++;
        I1 += ncal;
    }
    *neach = ncal;
    I0 = I1 - ncal;
    I1--;
    *Iinit = I0;
    *Ifin = I1;
}

//----------------- read elements from dump file --------------------//
/* Output: velo[Nk*Nrot*Nmd] that includes only data for SWNT
 *
 * SEDTYPE:
 * 0.x, 1.y, 2.z, 3.r, 4.theta
 */
void get_ROT_INFO(int ia, int nuat, int Nrot, int *irot, int *isite){
    *irot  = (ia%nuat)%Nrot;
    *isite = int((ia%nuat)/Nrot);
}
void get_ROT_INFO2(int iatom, int nuat, int Nrot, int *iz, int *irot, int *isite){
    *iz = int(iatom/nuat);
    *irot = (iatom%nuat)%Nrot;
    *isite = int((iatom%nuat)/Nrot);
}
int get_irow2read(const std::string &str1, const std::string &sword){
    int n = 0;
    std::string str2;
    std::stringstream ss(str1);
    while(std::getline(ss, str2, ' ')){
        if(str2 == "\0") continue;
        if(str2 == sword) return n - 2;
        n++;
    }
    std::cout << "Error: cannot find " << sword << std::endl;
    return n;
}
double extract_data_from_line(const std::string &str1, int Iget){
    int n = 0;
    std::string str2;
    std::stringstream ss(str1);
    while(std::getline(ss, str2, ' ')){
        if(str2 == "\0") continue;
        if(n==Iget) break;
        n++;
    }
    //std::cout << str1 << ", " << Iget << ", "  << str2 << std::endl;
    return std::stof(str2);
}

/**
 * Extract atomic velocities for a given atom and a given direction from dump file 
 *
 * Parameters
 * ----------
 *  filename : dump file name
 *  nmd : total number of MD steps
 *  natoms : number of atoms
 *  atom_index : atomic index, default=0
 *  label2read : label in dump file such as vx, vy, or vz
 *  nskip : number of columns to be skipped
 *
 * Return
 * ---------
 *  velocities : atomis velocities, shape=(nmd)
 *
 */
void extract_velocities_each_atom(
        char* filename, double* velocities, 
        int nmd, int natoms,
        int atom_index, char* label2read,
        const int nskip
        )
{
    int i, j, ia, count;
    int idx_md;
    char line[1024];
    double values[3];
    
    // --- get row index to read
    std::ifstream ifs0(filename);
    for(i=0; i<nskip; i++)
        ifs0.getline(line, sizeof(line));
    int irow2read = get_irow2read(line, label2read);
    ifs0.close();
    
    // --- read velocities
    std::ifstream ifs(filename);
    count = 0;
    for( idx_md=0; idx_md<nmd; idx_md++ ){

        //--- skip initial lines
        for(i=0; i<nskip; i++)
            ifs.getline(line, sizeof(line));
        
        //--- for all atom at a block
        for(ia=0;ia<natoms;ia++){
            
            ifs.getline(line, sizeof(line));
            
            // --- read data of ``ia`` atom at ``idx_md`` step.
            if ( ia == atom_index ){
                //cout << idx_md << " " << line << endl;
                velocities[idx_md] = extract_data_from_line(line, irow2read);
            }
        }
    }
    ifs.close();
}


//---- sed[nk*nc*nw]
void output_sed(char *krotfile, char *sumfile, char *LABEL, int nk, int nc, int nw, double *wlist, double *sed){
    int ik, ic, iw, num;
    double kz, kc;
    FILE *fp1 = fopen(krotfile, "w");
    FILE *fp2 = fopen(sumfile, "w");
    if(strcmp(LABEL, "T")==0) sprintf(LABEL, "theta");
    fprintf(fp1, "#Wavevector(k) Wavevector(theta) Frequency[THz]    SED(%s)\n", LABEL);
    fprintf(fp2, "#Wavevector(k)  Frequency[THz]   SED(%s)\n", LABEL);
    
    //--- for SED(kz,w)
    double *sedw = new double[nw];
    
    for(ik=0; ik<nk; ik++){
        kz = double(ik) / double(nk - 1.);
        for(iw=0; iw<nw; iw++) sedw[iw] = 0.;
        for(ic=0; ic<nc; ic++){
            kc = double(ic) / double(nc - 1.);
            for(iw=0; iw<nw; iw++){
                //--- for SED(kz, kt, w)
                num = ik*nc*nw + ic*nw + iw;
                fprintf(fp1, "%13.7f %13.7f   %15.10f %15.7e\n", kz, kc, wlist[iw], sed[num]);
                sedw[iw] += sed[num];
            }
            if(ik!=nk-1 || ic!=nc-1) fprintf(fp1, "\n");
        }
        //--- for SED(kz, w)
        for(iw=0; iw<nw; iw++){
            fprintf(fp2, "%13.7f %15.10f %15.7e\n", kz, wlist[iw], sedw[iw]);
        }
        if(ik!=nk-1) fprintf(fp2, "\n");
    }
    printf("Output: %s and %s\n", krotfile, sumfile);
    fclose(fp1);
    fclose(fp2);
    delete[] sedw; sedw = NULL;
}

void mkinput_sed(char *file){
    FILE *fp = fopen(file, "w");
    fprintf(fp, "DUMP_FILE   md4band.dump\n");
    fprintf(fp, "INIT_COORD  data.lammps\n");
    fprintf(fp, "NATOM_UNIT  40\n");
    fprintf(fp, "N_ROTATION  10\n");
    fprintf(fp, "\n");
    fprintf(fp, "TIMESTEP    0.0005    #-- [ps]\n");

    fclose(fp);
}

void conv_xyz2rtz(int nat, double **xyz0, double **rtz0, double *center){
    int i, j, ia;
    double length, theta;
    int Icoord[3];
    for(j=0;j<3;j++) Icoord[j] = j;
    double const pi = 3.14159265359;

    double **xyz = new double *[nat];
    for(i=0;i<nat;i++) xyz[i] = new double[3];
    for(i=0;i<nat;i++) for(j=0;j<3;j++) xyz[i][j] = xyz0[i][j] - center[j];
    
    //--- z-coord
    for(i=0;i<nat;i++) rtz0[i][2] = xyz[i][Icoord[2]];

    //--- r-coord
    for(ia=0; ia<nat; ia++){
        length = sqrt(xyz[ia][Icoord[0]]*xyz[ia][Icoord[0]] + xyz[ia][Icoord[1]]*xyz[ia][Icoord[1]]);
        rtz0[ia][0] = length;
    }

    //--- theta-coord
    for(ia=0; ia<nat; ia++){
        if(fabs(xyz[ia][Icoord[0]]) < 1e-10){
            if(xyz[ia][Icoord[1]] > 0.) theta = pi * 0.5;       
            else if(xyz[ia][Icoord[1]] <= 0.) theta = pi * 1.5;       
        }
        else{
            theta = atan(xyz[ia][Icoord[0]] / xyz[ia][Icoord[1]]);
            if(theta < 0.) theta = 2.*pi + theta;
            if(xyz[ia][Icoord[0]] < 0.){
                if(xyz[ia][Icoord[1]] > 0.) theta -= pi;
                else theta += pi;
            }
        }
        rtz0[ia][1] = theta;
    }

    //--- delte matrix 
    for(i=0;i<nat;i++){
        delete[] xyz[i]; xyz[i] = NULL;
    }
    delete[] xyz; xyz = NULL;
}

void substract_cylind(int Nat, double **vrtz, double **rtz0, double **rtz1){
    int i, j;
    const double pi = 3.14159265359;
    for(i=0;i<Nat;i++){
        for(j=0;j<3;j++){
            vrtz[i][j] = rtz1[i][j] - rtz0[i][j];
        }
        if(fabs(vrtz[i][1]) > pi) vrtz[i][1] = (vrtz[i][1]/fabs(vrtz[i][1])) * (- 2.*pi + vrtz[i][1]);
    }
}

