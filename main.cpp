#include <QCoreApplication>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <cmath>
#include <math.h>
using namespace std;

void SHA31(int darr[1600], int dddarr[5][5][64]);       //functions below
void SHA32(int dddarr[5][5][64], int darr[1600]);
void theta(int ain[5][5][64], int aout[5][5][64]);
void rho(int ain[5][5][64], int aout[5][5][64]);
void pi(int ain[5][5][64], int aout[5][5][64]);
void chi(int ain[5][5][64], int aout[5][5][64]);
void iota(int ain[5][5][64], int aout[5][5][64], int round);

int main()                          //implements SHA3
{
    int v[1600];                    //string to perform SHA3 on
    v[0] = 1;                       //gives both the 1 and 0 padding
    v[1087] = 1;
    for (int i = 1; i <1087;i++){
        v[i] = 0;
    }
    for (int i = 1088; i < 1600;i++){
        v[i] = 0;
    }
    int aoutddd1[5][5][64];         //3d arrays to hold the inputs and outputs
    int aoutddd2[5][5][64];         //updated by each function

    SHA31(v, aoutddd1);             //converts from 1 d array to 3d array
    for (int j = 0; j < 24; j++){
        theta(aoutddd1,aoutddd2);   //runs through each function, taking the previous output as input
        rho(aoutddd2,aoutddd1);
        pi(aoutddd1,aoutddd2);
        chi(aoutddd2,aoutddd1);
        iota(aoutddd1,aoutddd2,j);

        j++;                        //iterate within to prevent ending on the second array
        theta(aoutddd2,aoutddd1);   //since the loop starts with the first
        rho(aoutddd1,aoutddd2);     //to prevent missing the last output
        pi(aoutddd2,aoutddd1);
        chi(aoutddd1,aoutddd2);
        iota(aoutddd2,aoutddd1,j);
    }
    int darr[1600];                 //to hold the entries in a 1d array
    SHA32(aoutddd1,darr);
    for (int i = 0; i < 256;i++){
        cout<<darr[i]<<" ";
    }
    return 0;
}


void SHA31(int darr[1600], int dddarr[5][5][64]){ //turns 1 dimensional array in 3 dimensional array
    for (int i = 0; i <= 4; i++){
        for (int j = 0; j <= 4; j++){
            for (int k= 0; k <= 63; k++){
                dddarr[i][j][k] = darr[64*(5*j+i)+k];
            }
        }
    }
}

void SHA32(int dddarr[5][5][64], int darr[1600]){ //turns 3 dimensional array into 1 dimensional array
    for (int i = 0; i <= 4; i++){
        for (int j = 0; j <= 4; j++){
            for (int k= 0; k <= 63; k++){
                darr[64*(5*j+i)+k]=dddarr[i][j][k];
            }
        }
    }
}

void theta(int ain[5][5][64], int aout[5][5][64]){ //applies θ
    for (int i = 0; i <= 4; i++){
        for (int j = 0; j <= 4; j++){
            for (int k= 0; k <= 63; k++){
                int y = 0,                  //y and z are the summations for the equation
                    z = 0;
                for (int jp = 0; jp<= 4; jp++){ //gives ∑4 j ′=0 ain[i − 1][j′][k] and ∑4 j ′=0 ain[i + 1][j′][k − 1]
                    y += ain[(i+4)%5][jp][k]; //i-1 is the same as (i+4)mod5
                    y=y%2;

                    z += ain[(i+1)%5][jp][(k+63)%64]; //k-1 is the same as (k+63)mod64
                    z=z%2;

                }
                aout[i][j][k] = (ain[i][j][k]+y+z)%2; //applies the XOR operation
            }
        }
    }
}

void rho(int ain[5][5][64], int aout[5][5][64]){    //shifts ij drawer ina  cyclical shift
    int rhomatrix[5][5] = {{0,36,3,41,18},          //(((t+1)*(t+2))/2)%64 can be found in row i and column j
                           {1,44,10,45,2},          //of the given matrix
                           {62,6,43,15,61},         //eliminating the need for t
                           {28,55,25,21,56},
                           {27,20,39,8,14}};
    for (int i = 0; i <=4; i++){
        for (int j = 0; j <=4; j++){
            for (int k= 0; k <=63; k++){
            aout[i][j][k] = ain[i][j][(k-rhomatrix[i][j]+64)%64];
            }
        }
    }
}

void pi(int ain[5][5][64], int aout[5][5][64]){     //permutes within each drawer/slice with fixed k
    for (int ip = 0; ip <=4; ip++){
        for (int jp = 0; jp <=4; jp++){
            for (int k= 0; k <=63; k++){
                aout[jp][(2*ip + 3*jp)%5][k] = ain[ip][jp][k];
            }
        }
    }
}

void chi(int ain[5][5][64], int aout[5][5][64]){    //provides nonlinearity for each slice
    for (int i = 0; i <=4; i++){
        for (int j = 0; j <=4; j++){
            for (int k= 0; k <=63; k++){            //addition and mod 2 fulfill XOR
                aout[i][j][k] = (ain[i][j][k] + (((ain[(i + 1)%5][j][k] + 1)%2)*(ain[(i + 2)%5][j][k])) )%2;
            }
        }
    }
}

void iota(int ain[5][5][64], int aout[5][5][64],int round){
    int rc[168] = {1, 0, 0, 0, 0, 0, 0, 0};                         //constant bits for x^t reduced
    for (int i = 0; i <=4; i++){                                    //in F2[x]/(x 8 + x 6 + x 5 + x 4 + 1)
        for (int j = 0; j <=4; j++){
            for (int k= 0; k <=63; k++){
                aout[i][j][k] = ain[i][j][k];                       //copies over array for rest of iota
            }
        }
    }
    for (int next = 8 ; next <168; next++){                         //use LSR for computing rc for all 24 rounds
        rc[next] = (rc[next-8]+rc[next-4]+rc[next-3]+rc[next-2])%2;
    }
    for (int l = 0; l <7; l++){
        aout[0][0][(int)pow(2,l) - 1] = (ain[0][0][(int)pow(2,l) - 1]+rc[l + 7*round])%2;
    }
}

