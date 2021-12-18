//
//  file_handling.h
//
//
//  Created by masato on 5/26/15.
//
//

#ifndef ____file_handling__
#define ____file_handling__

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string.h>
using namespace std;

namespace file{
    void searchword(char *file, char *sword, char *data);
    int getNword(char line[]);
    void file_check(char file[]);
    int get_num_of_lines(char *file);
}

#endif /* defined(____file_handling__) */
