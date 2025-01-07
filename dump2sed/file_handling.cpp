//
//  file_handling.cpp
//
//  Created by masato on 5/26/15.
//
#include "file_handling.h"

namespace file{
  //
  // file: file name
  // sword: searched word
  // data: data which can be found just after "sword"
  void searchword(char *file, char *sword, char *data){
    int i, j, Ndata, check = 0;
    char line[1024], data0[128], *temp;
    ifstream ifs(file);
    while(!ifs.eof()){
      ifs.getline(line, sizeof(line));
      if(line[0]!='#'){
        Ndata = getNword(line);
        if(Ndata>1){
          temp = strtok(line, " "); sprintf(data0, temp);
          if(strcasecmp(data0, sword)==0){
            temp = strtok( NULL, " "); sprintf(data, temp);
            check = 1; 
            break;
          }
        }
      }
    }
    if(check==0){
        printf("cannot find %s in %s or no data\n", sword, file); 
        exit(0);
    }
  }
  //
  int getNword(char line[]){
    int Nw = 0; char *temp, line2[128];
    sprintf(line2, line); temp = strtok(line2, " ");
    while(temp!=NULL){ temp = strtok( NULL, " "); Nw++; }
    return Nw;
  }
  //
  void file_check(char file[]){
    ifstream ifs(file);
    if(!ifs){
      printf(" No such file exit. %s\n", file);
      exit(0);
    }
  }

  int get_num_of_lines(char *file){
    int count = 0;
    ifstream ifs(file);
    if(ifs){
        string line;
        while(true){
            getline(ifs, line);
            count++;
            if(ifs.eof()) break;
        }
    }
    ifs.close();
    return count - 1;
  }
}
