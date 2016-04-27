
#include "readinputsvl.h"
#include <iostream>
#include <fstream>
#include <string.h>

using namespace std;

// read the critical points information from the input data
void readinputsvl(char *infile,
                  svVector3Array *vec3profile,
                  svVector3Array *vec3oriProfile,
                  int *seed_num)
{
    
    if(infile == NULL) {
        //cerr << "inputfile doesn't exit .. do nothing..." << endl;
        return;
    };
    
    // read from file to storage
    ifstream inf(infile);
    char lineptr[256];
    double x, y, z, cx,cy,cz,ca;
    int counter, i;
    svInt  lineIndex=0;
    svScalar    lw;
    svVector3   ori;
    svVector3   pos;
    
    // get the seed number
    inf.getline(lineptr,  256);
    sscanf(lineptr, "%d", &counter);
    *seed_num = counter;
    
    while(1)
    {
        if(inf.getline(lineptr,  256))
        {
            int cc = 0;
            sscanf(lineptr, "%d", &counter);
            if(counter==1) {
                (*seed_num)--;
                inf.getline(lineptr,  256);
                continue;
            }
            else
            {
                for(i=0; i<counter; i++)
                {
                    inf.getline(lineptr, 256);
                    sscanf(lineptr, "%lf %lf %lf %lf %lf %lf",
                           &x, &y, &z,
                           &cx, &cy, &cz);
                    {
                        pos = svVector3(x,y,z);
                        ori = svVector3(cx,cy,cz);
                        
                        vec3profile[lineIndex].add(pos);
                        vec3oriProfile[lineIndex].add(ori);
                    }; // ++cc
                }; // loop of counter
                lineIndex++;
            }; // end if(counter)
        } else break; // end if(inf.getline)
    }; // end while(1)
    inf.close();
}

