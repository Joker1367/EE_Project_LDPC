#include <iostream>
#include <stdlib.h>
#include <cmath>

using namespace std;

const long long int para_1 = 4101842887655102017LL;
const long long int para_2 = 2685821657736338717LL;
const double para_3 = 5.42101086242752217E-20;
double STD_DEV;

unsigned long long int SEED = 3;
unsigned long long int RANV; 
int RANI = 0;

double Ranq1(){
    if(RANI == 0){
        RANV = SEED ^ para_1;
        RANV ^= RANV >> 21;
        RANV ^= RANV << 35;
        RANV ^= RANV >> 4;
        RANV *= para_2;
        RANI++;
    }
    RANV ^= RANV >> 21;
    RANV ^= RANV << 35;
    RANV ^= RANV >> 4;
    return RANV * para_2 * para_3;
}

void Normal(double& n1, double& n2, double std_dev){
    double x1, x2, s;
    do{
        x1 = Ranq1();
        x2 = Ranq1();
        x1 = 2 * x1 - 1;
        x2 = 2 * x2 - 1;
        s = x1 * x1 + x2 * x2;
    } while ( s >= 1.0);
    n1 = std_dev * x1 * sqrt(-2 * log(s) / s);
    n2 = std_dev * x2 * sqrt(-2 * log(s) / s);
}
