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

void extract_elements_2D_double(double **Morig, int Nget, int Ndat, int *Iget, double **Mget);
void cal_Irange_MPI(int *Iinit, int *Ifin, int *neach, int myrank, int nprocs, int Nloop);

void mkinput_sed(char *file);

//--- for read dump file
void get_ROT_INFO(int ia, int nuat, int Nrot, int *irot, int *isite);
void get_ROT_INFO2(int iatom, int nuat, int Nrot, int *iz, int *irot, int *isite);
int get_iline2read(const std::string &str1, const std::string &sword);
double extract_data_from_line(const std::string &str1, int Iget);
void extract_isite_velo_for_fft3D(char *FDUMP, int ISITE_GET, int SEDTYPE, int ITYPE_CNT, 
        int Nall, int Nz, int Nrot, int Nmd, int Nsite, double *velo_site, double radius);
void conv_xyz2rtz(int nat, double **xyz0, double **rtz0, double *center);
void substract_cylind(int Nat, double **vrtz, double **rtz0, double **rtz1);

//--- results
void output_sed(char *krotfile, char *sumfile, char *LABEL, int nk, int nc, int nw, double *wlist, double *sed);

int main(int argc, char **argv){
    
    //---- start MPI
    int nprocs, myrank;
    int Iinit, Ifin, neach;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if(myrank==0){
        printf("---------------------------------\n");
        printf("\n");
        printf("       dump2sed ver.1.0\n");
        printf("\n");
        printf("---------------------------------\n");
    }
    
    //--- parameters
    int i, j, ii, ia, count;
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
    file::searchword(IFILE, "DIRECTION", line); ised = atoi(line);
    file::searchword(IFILE, "DUMP_FILE", FDUMP);
    file::searchword(IFILE, "INIT_COORD", FXYZ);
    file::searchword(IFILE, "NATOM_UNIT", line); Nuat = atoi(line);
    file::searchword(IFILE, "N_ROTATION", line); Nrot = atoi(line);
    file::searchword(IFILE, "TIMESTEP", line);   tstep = atof(line);
    int ITYPE_CNT = 1;

    //--- read initial xyz file
    int Natom = lmp::get_natom_from_input(FXYZ);
    int Ntype = lmp::get_ntype_from_input(FXYZ);
    int *ID_all = new int[Natom];        //--------- new ID_ALL[Natom]
    int *TY_all = new int[Natom];        //--------- new TY_all[Natom]
    double cube[6];
    double *mass = new double[Ntype];    //--------- new mass[Ntype]
    double **Call = new double *[Natom];
    for(i=0;i<Natom;i++) Call[i] = new double[3];
    lmp::get_coord_from_input(FXYZ, Ntype, mass, Natom, ID_all, TY_all, cube, Call);
    
    //--- extract tube data: 
    int ID_SWNT_at_MOL = 1;
    int Ncnt = 0;
    for(i=0;i<Natom;i++){
        if(TY_all[i] == ID_SWNT_at_MOL) Ncnt++;
    }
    int *I_cnt = new int[Ncnt];            //---- Ccnt[i][j] == Call[I_cnt[i]][j]
    double **Ccnt = new double *[Ncnt];    //---- Ccnt[Ncnt][3]
    for(i=0;i<Ncnt;i++) Ccnt[i] = new double[3];
    count = 0;
    for(ia=0;ia<Natom;ia++){
        if(TY_all[ia] == ID_SWNT_at_MOL){
            I_cnt[count] = ia;
            count++;
        }
    }
    extract_elements_2D_double(Call, Ncnt, 3, I_cnt, Ccnt); //--- Ccnt = Call[I_cnt]

    //--- radius
    int XIND = 0, YIND = 1;
    double radius = 0., center[3], xx, yy;
    for(i=0;i<3;i++) center[i] = 0.;
    for(i=0;i<Ncnt;i++){
        for(j=0;j<3;j++) center[j] += Ccnt[i][j] / double(Ncnt);
    }
    for(i=0;i<Ncnt;i++){
        xx = Ccnt[i][XIND] - center[0];
        yy = Ccnt[i][YIND] - center[1];
        radius += sqrt(xx*xx + yy*yy) / double(Ncnt);
    }

    //--- delete Call
    for(i=0;i<Natom;i++){
        delete[] Call[i];
        Call[i] = 0;
    }
    delete[] Call; Call = 0;
    delete[] ID_all; ID_all = 0;
    delete[] TY_all; TY_all = 0;

    //--- read dump file
    int nstep, Nmd, nline = 0;
    nstep = lmp::get_dumpstep(FDUMP);
    nline = file::get_num_of_lines(FDUMP);
    Nmd   = nline / (Natom + 9);
    Nmd   = pow(2, int(log2(Nmd)));    //--- Nmd = 2^**
    
    //--- Frequency [THz]
    double interval = tstep * nstep;                   //--- dumping interval [ps]
    double ttotal   = double(Nmd) * double(interval);  //--- total duration [ps]
    double wmax     = 1. / interval;
    double *wlist = new double[Nmd];
    for(i=0;i<Nmd;i++) wlist[i] = wmax * double(i) / double(Nmd-1);
    
    //--- for parameters after FFT
    int Mcnt = int(Ncnt/Nuat);   //--- # of unit cells
    int Nsite = int(Nuat/Nrot);  //--- # of sites
    int nk = int(Mcnt/2) + 1;        //--- # of k_z
    int nc = int(Nrot/2) + 1;        //--- # of k_theta
    int nw = int(Nmd/2) + 1;         //--- # of frequency

    //--- output variables
    if(myrank == 0){
        printf("============= variables ==============\n");
        printf(" Step time [ps]       : %7.5f\n", tstep);
        printf(" Interval  [ps]       : %7.5f\n", interval);
        printf(" Max. freq [THz]      : %5.2f\n", wmax);
        printf(" # of MD steps        : %3d\n", Nmd);
        printf(" # of all atoms       : %3d\n", Natom);
        printf(" # of atoms at SWNT   : %3d\n", Ncnt);
        printf("   = %d atoms/unit x %d units\n", Nuat, Mcnt);
        printf("             (kz:%2d, kt:%2d, w:%3d)\n", nk, nc, nw);
        printf(" Radius [A]           : %7.4f\n", radius);
        printf("======================================\n");
        printf("\n");
    }

    //----- velo[Nrot*Nk*Nmd]
    /*
     * SEDTYPE:
     * 0. x, 1.y, 2.z, 3.r, 4.theta
     *
     */
    int N_sedtype = 5;
    
    cal_Irange_MPI(&Iinit, &Ifin, &neach, myrank, nprocs, N_sedtype); 
    if(ised >= 0 and ised <= 4){
        Iinit = ised;
        Ifin = ised;
    }
    
    int SEDTYPE, isite, ik, ic, iw, ik2, ic2, iw2;
    double spectral;
    double *velo = new double[Mcnt*Nrot*Nmd];           //---- Large array: velo[Mcnt*Nrot*Nmd]
    double *sed_raw = new double[Mcnt*Nrot*Nmd];        //---- Large array: sed_raw[Mcnt*Nrot*Nmd]
    double *sed = new double[nk*nc*nw];
    char OFILE1[256], OFILE2[256], LABEL[128];
    for(i=0; i<nk*nc*nw; i++) sed[i] = 0.;
    
    //printf("myrand(%d): %d to %d\n", myrank, Iinit, Ifin);

    for(SEDTYPE = Iinit; SEDTYPE < Ifin+1; SEDTYPE++)
    {
        for(isite=0; isite<Nsite; isite++){

            if(myrank==0) printf("=== Start SITE %d/%d\n", isite+1, Nsite);
            
            //----- extract velocity and perform Fourier transform
            extract_isite_velo_for_fft3D(FDUMP, isite, SEDTYPE, ITYPE_CNT, 
                    Natom, Mcnt, Nrot, Nmd, Nsite, velo, radius);
            MYFFTW::fftw3d(Mcnt, Nrot, Nmd, velo, sed_raw);
            for(i=0;i<Mcnt*Nrot*Nmd;i++) sed_raw[i] /= double(Mcnt * Nrot * Nmd);
            
            for(ik=0; ik<nk; ik++){
                for(ic=0; ic<nc; ic++){
                    for(iw=0; iw<nw; iw++){
                        spectral = sed_raw[ik*Nrot*Nmd + ic*Nmd + iw];
                        sed[ik*nc*nw + ic*nw + iw] += spectral * 0.5;
                    }
                }
            }
        }
        if(SEDTYPE == 0){
            sprintf(LABEL, "X");
        }
        else if(SEDTYPE == 1){
            sprintf(LABEL, "Y");
        }
        else if(SEDTYPE == 2){
            sprintf(LABEL, "Z");
        }
        else if(SEDTYPE == 3){
            sprintf(LABEL, "R");
        }
        else if(SEDTYPE == 4){
            sprintf(LABEL, "T");
        }
        else{
            printf("Error: SET TYPE %d is not defined.\n", SEDTYPE);
            exit(0);
        }
        
        //---- output SED file!
        sprintf(OFILE1, "SED_kt_%s.txt", LABEL);
        sprintf(OFILE2, "SED_%s.txt", LABEL);
        output_sed(OFILE1, OFILE2, LABEL, nk, nc, nw, wlist, sed);
    }

    delete[] velo;
    MPI_Finalize();
    return 0;
}

//--- extract Iget[Nget] elements from Morig[**][Ndat] and output Mget[Nget][Ndat]
void extract_elements_2D_double(double **Morig, int Nget, int Ndat, int *Iget, double **Mget){
    int i, j;
    for(i=0; i<Nget; i++){
        for(j=0; j<Ndat; j++){
            Mget[i][j] = Morig[Iget[i]][j];
        }
    }    
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
int get_iline2read(const std::string &str1, const std::string &sword){
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
void extract_isite_velo_for_fft3D(char *FDUMP, int ISITE_GET, int SEDTYPE, int ITYPE_CNT,
        int Nall, int Nz, int Nrot, int Nmd, int Nsite, double *velo_site, double radius)
{
    int i, j, ir, ia, iua, num, count, count_cnt, count_read, itype;
    int I_MD, I_Z, I_ROT, I_SITE;
    int ILINE = -1;
    int nuat = Nsite * Nrot;
    char line[1024], LABEL[16];
    int Icoord[3], Ivelo[3];
    double center[3], coord[3];

    const int Natom = Nsite * Nrot * Nz;
    const int Nskip = 9;

    //------ get line numbers to read
    std::ifstream ifs0(FDUMP);
    for(i=0;i<Nskip;i++) ifs0.getline(line, sizeof(line));
    Icoord[0] = get_iline2read(line, "x");
    Icoord[1] = get_iline2read(line, "y");
    Icoord[2] = get_iline2read(line, "z");
    Ivelo[0]  = get_iline2read(line, "vx");
    Ivelo[1]  = get_iline2read(line, "vy");
    Ivelo[2]  = get_iline2read(line, "vz");
    ifs0.close();
    
    //
    //--- prepare for reading velocity (info. will be read for each unit.)
    //
    double **xyz_unit = new double *[Nrot];            //--- it's for the conversoin to cylindrical coord.
    for(i=0;i<Nrot;i++) xyz_unit[i]  = new double[3];
    double **vxyz_unit = new double *[Nrot];           //--- velocity
    for(i=0;i<Nrot;i++) vxyz_unit[i] = new double[3];
    double **rtz0_unit = new double *[Nrot];           //--- rtz: before (t = t0)
    for(i=0;i<Nrot;i++) rtz0_unit[i] = new double[3];
    double **rtz1_unit = new double *[Nrot];           //--- rtz: after  (t = t0 + dt)
    for(i=0;i<Nrot;i++) rtz1_unit[i] = new double[3];
    double **vrtz_unit = new double *[Nrot];           //--- vrtz: after  (t = t0 + dt)
    for(i=0;i<Nrot;i++) vrtz_unit[i] = new double[3];
    
    //
    //--- read velocities!! velo_site[Nz*Nrot*Nmd]
    //
    for(i=0;i<Nz*Nrot*Nmd;i++) velo_site[i] = 0.;

    std::ifstream ifs(FDUMP); count = 0;
    for( I_MD=0; I_MD<Nmd; I_MD++ ){

        //--- skip initial lines
        for(i=0; i<Nskip; i++){
            ifs.getline(line, sizeof(line));
        }
        
        //--- for all atom at a block
        count_cnt = 0;
        count_read = 0;
        for(ia=0;ia<Nall;ia++){
            ifs.getline(line, sizeof(line));
            itype = int(extract_data_from_line(line, 1));
            if (itype != ITYPE_CNT) continue;       //--- C60 info. is skipped.

            get_ROT_INFO2(count_cnt, nuat, Nrot, &I_Z, &I_ROT, &I_SITE);
            
            //--- center of unit cell
            for(j=0;j<3;j++) coord[j] = extract_data_from_line(line, Icoord[j]);
            for(j=0;j<3;j++) center[j] += coord[j] / double(nuat);
             
            //--- read only data for {ISITE_GET} atom
            if( I_SITE == ISITE_GET ){
                for(j=0;j<3;j++){
                    xyz_unit[I_ROT][j] = coord[j];
                    vxyz_unit[I_ROT][j] = extract_data_from_line(line, Ivelo[j]);
                }
                count_read++;
            }
            
            count_cnt++;
            
            //--- extract velocities for each unit cell
            if( count_cnt % nuat == 0 ){
                if(count_read%Nrot != 0){
                    printf("Wired!!!!!!\n");
                    exit(0);
                }
                
                //
                //--- 2. convert {xyz_unit, vxyz_unit} to {velo_site} 
                //
                //--- 2a. for (x,y,z)
                //
                if(SEDTYPE < 3){ 
                    for(ir=0; ir<Nrot; ir++){
                        num = I_Z*Nrot*Nmd + ir*Nmd + I_MD;
                        if(num>=Nz*Nrot*Nmd){
                            printf("Error too large\n");
                            exit(0);
                        }
                        velo_site[num] = vxyz_unit[ir][SEDTYPE];
                        count++;
                    }
                }
                //
                //--- 2b. for (r,theta,z). Conversion to cylindrical coord is needed.
                //--- Input:  xyz_unit[Nrot][3], vxyz_unit[Nrot][3], center[3]
                //--- ====>:  rtz0_unit[Nrot][3], rtz1_unit[Nrot][3]
                //--- Output: vrtz_unit[Nrot][3]
                //
                else if(SEDTYPE < 5){
                    //--- rtz(t0)
                    conv_xyz2rtz(Nrot,  xyz_unit, rtz0_unit, center);
                    //--- xyz(t0 + dt)
                    for(i=0;i<Nrot;i++){
                        for(j=0;j<3;j++) xyz_unit[i][j] += vxyz_unit[i][j];
                    }
                    //--- rtz(t0 + dt)
                    conv_xyz2rtz(Nrot, xyz_unit, rtz1_unit, center);
                    //--- vrtz(t0)
                    substract_cylind(Nrot, vrtz_unit, rtz0_unit, rtz1_unit);
                    for(i=0;i<Nrot;i++) vrtz_unit[i][1] *= radius;
                    
                    for(ir=0; ir<Nrot; ir++){
                        num = I_Z*Nrot*Nmd + ir*Nmd + I_MD;
                        velo_site[num] = vrtz_unit[ir][SEDTYPE-3];
                        count++;
                    }
                }
                else{
                    printf("Error");
                    exit(0);
                }
                
                //--- initialize
                for(j=0;j<3;j++) center[j] = 0.;
            }    
        }
    } 
    ifs.close();
    if(count != Nz*Nrot*Nmd ){
        printf("Error: %d != %d\n", count, Nz*Nrot*Nmd);
        exit(0);
    }

    for(i=0;i<Nrot;i++) delete[] xyz_unit[i];
    delete[] xyz_unit; xyz_unit = NULL;
    for(i=0;i<Nrot;i++) delete[] vxyz_unit[i];
    delete[] vxyz_unit; vxyz_unit = NULL;
    for(i=0;i<Nrot;i++) delete[] rtz0_unit[i];
    delete[] rtz0_unit; rtz0_unit = NULL;
    for(i=0;i<Nrot;i++) delete[] rtz1_unit[i];
    delete[] rtz1_unit; rtz1_unit = NULL;
    for(i=0;i<Nrot;i++) delete[] vrtz_unit[i];
    delete[] vrtz_unit; vrtz_unit = NULL;
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

