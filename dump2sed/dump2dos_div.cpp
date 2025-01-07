#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <math.h>
#include <iostream>
#include <omp.h>
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

void mkinput_dos(char *file);

// --- to read dump file
int get_irow2read(const std::string &str1, const std::string &sword);
double extract_data_from_line(const std::string &str1, int Iget);

void extract_velocities_each_atom(
        char* filename, double* velocities,
        int nmd, int natoms,
        int atom_index = 0, char* label2read = "vx", const int nskip=9);

void print_array(int ndat, double* array);

// --- output the result
void write_dos_div(
        char* outfile, int nfreq, int ndiv, 
        double* frequencies, double** dos, double ldiv,
        int natoms, int* natoms_sec
        );

int main(int argc, char **argv){
    
    //if(myrank==0){
    printf("---------------------------------\n");
    printf("\n");
    printf("       dump2sed ver.1.0\n");
    printf("\n");
    printf("---------------------------------\n");
    //}
    
    // --- parameters
    int i, j, ii, ia, idat, count;
    int Nuat, Nrot, i_coord_type;
    char IFILE[128], word[128];
    char FDUMP[128], CTYPE[128], file_xyz[128];
    char line[256];
    
    // --- input file name
    if (argc > 1) sprintf(IFILE, "%s", argv[1]);
    else{
        sprintf(IFILE, "example_sed.in");
        printf("Input file name.\nSee %s\n", IFILE);
        mkinput_dos(IFILE);
        exit(0);
    }
    
    // --- file check
    ifstream ifs_check(IFILE);
    if(ifs_check.fail()){
        printf("Cannot find %s\n", IFILE);
        exit(0);
    }
    ifs_check.close();
    
    // --- read input file
    file::searchword(IFILE, "DUMP_FILE", FDUMP);
    //file::searchword(IFILE, "INIT_COORD", file_xyz);    
    file::searchword(IFILE, "TIMESTEP", line);
    double tstep = atof(line);
    file::searchword(IFILE, "NUM_DIVISIONS", line);
    int num_div = atoi(line);
    
    /**************
     * Read FDUMP
     ***************/
    // --- read initial xyz file
    int Natom = lmp::get_natom_from_dump(FDUMP);
    
    // --- read dump file
    int nstep = lmp::get_dumpstep(FDUMP);
    int nline = file::get_num_of_lines(FDUMP);
    int Nmd = nline / (Natom + 9);   // total number of output data
    Nmd = pow(2, int(log2(Nmd)));    // Nmd = 2^** data will be used
    
    // --- list of frequencies
    int nfreq;
    if(Nmd%2)
        nfreq = Nmd / 2;
    else
        nfreq = (Nmd + 1) / 2;
    
    double time_dump = tstep * nstep; // ps
    double wmax = 0.5 / time_dump;    // THz
    double delta_w = wmax / nfreq;      // THz
    double* frequencies = new double[nfreq];
    for(i=0; i<nfreq; i++) 
        frequencies[i] = delta_w * i;
    
    printf("\n");
    printf(" Number of data  : %d\n", Nmd);
    printf(" Number of atoms : %d\n", Natom);
    printf("\n");
    
    /*****************
     * set sections
     *****************/
    //// --- read cell size
    //int Ntype = lmp::get_ntype_from_input(file_xyz);
    //int *ID_all = new int[Natom];
    //int *TY_all = new int[Natom];
    int *idx_sec = new int[Natom];
    //double *mass = new double[Ntype];
    double cube[6];
    double **positions = new double *[Natom];
    for(i=0;i<Natom;i++)
        positions[i] = new double[3];
    
    //lmp::get_coord_from_input(file_xyz, Ntype, mass, Natom, ID_all, TY_all, cube, positions);
    lmp::read_init_structure_dump(FDUMP, cube, positions);

    // --- assign atoms to sections
    double ldiv = (cube[1] - cube[0]) / double(num_div);
    int *natoms_sec = new int[num_div];
    for(i=0; i<num_div; i++)
        natoms_sec[i] = 0;
    
    for(i=0; i<Natom; i++){
        j = int((positions[i][0] - cube[0]) / ldiv);
        if(j == num_div)
            j = num_div - 1;
        idx_sec[i] = j;
        natoms_sec[j] += 1;
    }
    
    /**********************
     * prepare arrays
     **********************/
    double *velo_each = new double[Nmd];
    double *dos_each = new double[Nmd];
    double **dos_total = new double *[num_div];
    for(i=0; i<num_div; i++)
        dos_total[i] = new double[Nmd];
    
    char labels[3][10] = {"vx", "vy", "vz"};
    
    for(i=0; i<num_div; i++)
        for(j=0; j<Nmd; j++)
            dos_total[i][j] = 0.;
    
    #pragma omp parallel for
    for(ia=0; ia<Natom; ia++){
        
        if(ia % 1000 == 0)
            printf(" %d/%d \n", ia, Natom);

        for(j=0; j<3; j++){
            
            extract_velocities_each_atom(FDUMP, velo_each, Nmd, Natom, ia, labels[j]);
            
            MYFFTW::fftw1d(Nmd, velo_each, dos_each);

            for(idat=0; idat<Nmd; idat++)
                dos_total[idx_sec[ia]][idat] += dos_each[idat] / double(Natom * 3.);
            
        }
    }
    // end of omp parallel
    
    // --- DOS
    int i1, i2;
    double** dos_ave = new double *[num_div];
    for(i=0; i<num_div; i++)
        dos_ave[i] = new double[nfreq];
    
    for (j=0; j<num_div; j++){
        for (i=0; i<nfreq; i++){
            i1 = i;
            i2 = Nmd - 1 - i;
            if (i1 > i2)
                break;
            double dos1 = dos_total[j][i1];
            double dos2 = dos_total[j][i2];
            dos_ave[j][i] = (dos1 + dos2) * 0.5;
        }
    }
    
    // --- output DOS
    char outfile_dos[100] = "dos_div.txt";
    write_dos_div(outfile_dos, nfreq, num_div, frequencies, dos_ave, ldiv, Natom, natoms_sec);
    
    delete[] velo_each;
    delete[] dos_each;
    delete[] frequencies;
    
    for(i=0; i<num_div; i++) delete[] dos_total[i];
    delete[] dos_total;
    dos_total = NULL;
    
    for(i=0; i<num_div; i++) delete[] dos_ave[i];
    delete[] dos_ave;
    dos_ave = NULL;
    
    return 0;
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

/***********************************************************************************
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

void mkinput_dos(char *file){
    FILE *fp = fopen(file, "w");
    fprintf(fp, "DUMP_FILE   md4band.dump\n");
    //fprintf(fp, "INIT_COORD  data.lammps\n");
    //fprintf(fp, "NATOM_UNIT  40\n");
    //fprintf(fp, "N_ROTATION  10\n");
    //fprintf(fp, "\n");
    fprintf(fp, "TIMESTEP    0.0005    #-- [ps]\n");
    fclose(fp);
}

void write_dos_div(
        char* outfile, int nfreq, int ndiv, 
        double* frequencies, double** dos, double ldiv,
        int natoms, int* natoms_sec)
{
    
    int idiv, iw;
    FILE* fp = fopen(outfile, "w");
    
    fprintf(fp, "# Number of sections : %d\n", ndiv);
    fprintf(fp, "# Length of sections : %.3f\n", ldiv);
    fprintf(fp, "# Number of atoms    : %d\n", natoms);
    fprintf(fp, "# Number of atoms in section\n");
    for(idiv=0; idiv<ndiv; idiv++){
        fprintf(fp, "#  section %d : %d \n", idiv+1, natoms_sec[idiv]);
    }
    fprintf(fp, "# Frequency[THz]   DOS(a.u.)1, 2, 3, ...\n");
    
    for(iw=0; iw<nfreq; iw++){
        fprintf(fp, "%15.10f  ", frequencies[iw]);
        for(idiv=0; idiv<ndiv; idiv++){
            fprintf(fp, "%13.8e ", dos[idiv][iw]);
        }
        fprintf(fp, "\n");
    }
    printf("Output: %s \n", outfile);
    fclose(fp);
}

void print_array(int ndat, double* array){
    for (int i=0; i<ndat; i++)
        printf(" %02d : %8.3f \n", i, array[i]);
    cout << endl;
}

