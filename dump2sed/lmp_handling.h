//
//  lmp_handling.h
//
//
//  Created by masato on 5/26/15.
//
//

#ifndef ____lmp_handling__
#define ____lmp_handling__

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string.h>

namespace lmp{

    //--- get # of types from data.lammps 
    int get_ntype_from_input(char *dfile){
        int i;
        char word[128];
        ifstream ifs(dfile);
        if(ifs.fail()){
            printf("Cannot find %s\n", dfile);
            exit(0);
        }
        for(i=0;i<4;i++) ifs >> word;
        ifs >> word;
        ifs.close();
        return atoi(word);
    }

    //--- get # of atoms from data.lammps
    int get_natom_from_input(char *dfile){
        char data[128];
        ifstream ifs(dfile);
        if(ifs.fail()){
            printf("Cannot find %s\n", dfile);
            exit(0);
        }
        ifs >> data >> data;
        ifs >> data;
        ifs.close();
        return atoi(data);
    }

    //--- get molecular info. from data.lammps
    //--- Input:  dfile, Ntype, Natom
    //--- Output: Ntype, mass, Natom, IDlist, TYPElist, cube, coord
    int get_coord_from_input(char *dfile, int Ntype, double *mass, 
            int Natom, int *ID_atom, int *TY_atom, double *cube, double **coord)
    {
        int i, j;
        char word[128];
        ifstream ifs(dfile);
        if(ifs.fail()){
            printf("Cannot find %s\n", dfile);
            exit(0);
        }
        for(i=0;i<7;i++) ifs >> word;
        if(strcmp(word, "types")!=0){
            printf("Error while reading %s (1)\n", dfile);
            exit(0);
        }
        
        //--- cube
        for(i=0;i<3;i++){
            ifs >> cube[2*i] >> cube[2*i+1];
            ifs >> word >> word;
        }
        
        //--- mass
        ifs >> word;
        for(i=0;i<Ntype;i++) ifs >> word >> mass[i] >> word;
        
        //--- atomic info.
        ifs >> word;
        if(strcmp(word, "Atoms")!=0){
            printf("Error while reading %s (2)\n", dfile);
            exit(0);
        }
        for(i=0;i<Natom;i++){
            ifs >> ID_atom[i] >> TY_atom[i];
            for(j=0;j<3;j++) ifs >> coord[i][j];
        }
        ifs.close();
        return 0;
    }
    
    //----------------- DUMP FILE -------------------//
    int get_dumpstep(char *fdump){
        int i, istep0, istep1, natom;
        char word[128], line[1024];
        ifstream ifs(fdump);
        if(ifs.fail()){
            printf("Cannot find %s\n", fdump);
            exit(0);
        }
        //--- first block
        ifs.getline(line, sizeof(line));
        ifs >> istep0;
        for(i=0;i<2;i++) ifs.getline(line, sizeof(line));
        ifs >> natom;
        ifs.getline(line, sizeof(line));
        for(i=0;i<natom+6;i++) ifs.getline(line, sizeof(line));
        ifs >> istep1;
        return istep1 - istep0;
    }
    
    //--- get # of atoms from dump file
    int get_natom_from_dump(char *filename){
        int i, istep0, istep1, natom;
        char word[128], line[1024];
        ifstream ifs(filename);
        if(ifs.fail()){
            printf("Cannot find %s\n", filename);
            exit(0);
        }
        //--- first block
        ifs.getline(line, sizeof(line));
        ifs >> istep0;
        for(i=0;i<2;i++)
            ifs.getline(line, sizeof(line));
        ifs >> natom;
        return natom;
    }

    //-- read initial structure in dump file
    int read_init_structure_dump(char *filename, double *cube, double **positions){
        
        int i, j;
        int idx, type;
        int timestep, natoms;
        double vx, vy, vz;
        char line[124];
        
        ifstream ifs(filename);
        if(ifs.fail()){
            printf("Cannot find %s\n", filename);
            exit(0);
        }
        // ITEM: TIMESTEP
        ifs.getline(line, sizeof(line));
        ifs >> timestep;
        
        // ITEM: NUMBER OF ATOMS
        ifs.getline(line, sizeof(line));
        ifs.getline(line, sizeof(line));
        ifs >> natoms;
        
        // ITEM: BOX BOUNDS pp pp pp
        ifs.getline(line, sizeof(line));
        ifs.getline(line, sizeof(line));
        for(j=0; j<3; j++)
            ifs >> cube[j*2] >> cube[j*2+1];
        
        // ITEM: ATOMS id type x y z vx vy vz
        ifs.getline(line, sizeof(line));
        ifs.getline(line, sizeof(line));
        for(i=0; i<natoms; i++){
            ifs >> idx >> type;
            for(j=0; j<3; j++)
                ifs >> positions[i][j];
            ifs >> vx >> vy >> vz;
        }

    }

}
#endif /* defined(____lmp_handling__) */

