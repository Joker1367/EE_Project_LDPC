#include <iostream>
#include <stdlib.h>
#include "RandVar.h"
#include "Graph.h"
#include <fstream>
#include <vector>
#include <string>

using namespace std;

vector<vector<int>> H;
vector<vector<int>> RRE_H;
vector<vector<int>> G;
vector<int> information;
vector<int> codeword;
vector<int> result;
vector<double> code;
vector<double> pre_sum;
vector<double> post_sum;

//int rotate_1[4][8] = {{-1, 440, 441, 420, -1, 359, 298, 196}, {413, -1, 38, 124, 162, -1, 448, 273}, {460, 21, -1, 41, 61, 102, -1, 265}, {48, 375, 423, -1, 299, 175, 13, -1}};
//int rotate_2[4][8] = {{-1, 229, 262, 30, -1, 322, 153, 14}, {260, -1, 332, 404, 275, -1, 32, 250}, {428, 232, -1, 431, 169, 139, -1, 447}, {201, 389, 129, -1, 186, 243, 429, -1}};
int rotate[4][8] = {{-1, 82, 170, 218, -1, 262, 150, 2}, {338, -1, 445, 4, 460, -1, 405, 292}, {250, 319, -1, 419, 38, 30, -1, 374}, {46, 419, 73, -1, 276, 220, 103, -1}};

int main(int argc, char* argv[]){
    int K, N, M, dc, dv, dc1, dc2;
    int e = 922;
    int position;
    int index;
    int threshold;
    bool correct = false;
    double bit_ep = 0;
    double R;
    string lines;
    ifstream inFile;
    ofstream outFile;

    // Read information of parity-check matrix
    outFile.open("database_QC_QAM.txt");
    N = 7376;
    M = 3688;
    dc = 3;
    dv = 6;

    Graph graph(N, M);

    // Read in Parity-Check Matrix
    H.assign(M, vector<int>(N, 0));
    RRE_H.assign(M, vector<int>(N, 0));

    // left-top
    for(int l = 0; l < 4; l++){
        for(int m = 0; m < 8; m++){
            int tmp = rotate[l][m];
            int row = l * e;
            int col = m * e;
            if(tmp == -1){
                for(int i = 0; i < e; i++){
                    for(int j = 0; j < e; j++){
                        H[i + row][j + col] = 0;
                        RRE_H[i + row][j + col] = 0;
                    }
                }
            }
            else{
                for(int i = 0; i < e; i++){
                    for(int j = 0; j < e; j++){
                        if(j == tmp){
                            H[i + row][j + col] = 1;
                            RRE_H[i + row][j + col] = 1;
                            graph.VAR[j + col].neighbor.push_back(i + row);
                            graph.CHK[i + row].neighbor.push_back(j + col);    
                        }
                        else{
                            H[i + row][j + col] = 0;
                            RRE_H[i + row][j + col] = 0;
                        }
                    }
                    tmp = (tmp + 1) % e;
                }
            }
        }
    }

    // right-top
    /*for(int l = 0; l < 4; l++){
        for(int m = 0; m < 8; m++){
            int tmp = rotate_2[l][m];
            int row = l * e;
            int col = m * e;
            if(tmp == -1){
                for(int i = 0; i < e; i++){
                    for(int j = 0; j < e; j++){
                        H[i + row][3688 + j + col] = 0;
                        RRE_H[i + row][3688 + j + col] = 0;
                    }
                }
            }
            else{
                for(int i = 0; i < e; i++){
                    for(int j = 0; j < e; j++){
                        if(j == tmp){
                            H[i + row][3688 + j + col] = 1;
                            RRE_H[i + row][3688 + j + col] = 1;
                            graph.VAR[3688 + j + col].neighbor.push_back(i + row);
                            graph.CHK[i + row].neighbor.push_back(3688 + j + col);    
                        }
                        else{
                            H[i + row][3688 + j + col] = 0;
                            RRE_H[i + row][3688 + j + col] = 0;
                        }
                    }
                    tmp = (tmp + 1) % e;
                }
            }
        }
    }

    // left-bottom
    for(int l = 0; l < 4; l++){
        for(int m = 0; m < 8; m++){
            int row = l * e;
            int col = m * e;
            for(int i = 0; i < e; i++){
                for(int j = 0; j < e; j++){
                    H[1844 + i + row][j + col] = 0;
                    RRE_H[1844 + i + row][j + col] = 0;
                }
            }
        }
    }

    // right-bottom
    for(int l = 0; l < 4; l++){
        for(int m = 0; m < 8; m++){
            int tmp = rotate_1[l][m];
            int row = l * e;
            int col = m * e;
            if(tmp == -1){
                for(int i = 0; i < e; i++){
                    for(int j = 0; j < e; j++){
                        H[1844 + i + row][3688 + j + col] = 0;
                        RRE_H[1844 + i + row][3688 + j + col] = 0;
                    }
                }
            }
            else{
                for(int i = 0; i < e; i++){
                    for(int j = 0; j < e; j++){
                        if(j == tmp){
                            H[1844 + i + row][3688 + j + col] = 1;
                            RRE_H[1844 + i + row][3688 + j + col] = 1;
                            graph.VAR[3688 + j + col].neighbor.push_back(1844 + i + row);
                            graph.CHK[1844 + i + row].neighbor.push_back(3688 + j + col);    
                        }
                        else{
                            H[1844 + i + row][3688 + j + col] = 0;
                            RRE_H[1844 + i + row][3688 + j + col] = 0;
                        }
                    }
                    tmp = (tmp + 1) % e;
                }
            }
        }
    }*/

////////////////////////////////////////////////////////////////////////////// Find Generator Matrix ///////////////////////////////////////////////////
    vector<int> change(N, -1);
    vector<int> pivot(M, -1);

    for(int i = 0; i < N; i++) change[i] = i;

    index = 0;

    // Gaussian elimination into row-echelon form
    for(int i = 0; i < N; i++){
        for(int j = index; j < M; j++){
            if(RRE_H[j][i] == 0) continue;
            else{
                for(int k = 0; k < N; k++) swap(RRE_H[index][k], RRE_H[j][k]);
                for(int k = index + 1; k < M; k++){
                    if(RRE_H[k][i] == 1) {
                        for(int l = 0; l < N; l++) RRE_H[k][l] ^= RRE_H[index][l];
                    }
                }
                pivot[index++] = i;
                break;
            }
        }
    }

    K = N - index;
    R = 1.0 * K / N;
    cout << R << endl;

    // initialize G
    G.assign(K, vector<int>(N, 0));

    cout << index << endl;

    // Reducued-row-ecelon form
    for(int i = 0; i < index; i++){
        for(int j = 0; j < i; j++){
            if(RRE_H[j][pivot[i]] == 1) for(int k = 0; k < N; k++) RRE_H[j][k] ^= RRE_H[i][k];
        }
    }

    // swap column of RRE_H into [A, I]
    int cur = N - 1;
    for(int i = index - 1; i >= 0; i--){
        int tmp = pivot[i];
        for(int j = 0; j < index; j++) swap(RRE_H[j][tmp], RRE_H[j][cur]);
        change[tmp] = change[cur];
        change[cur] = tmp;
        cur--;
    }

    // obtained G' by G' = [I, -A^T]
    for(int i = 0; i < N; i++){
        if(i < K){
            for(int j = 0; j < K; j++) G[j][i] = (j == i) ? 1 : 0;
        }
        else{
            for(int j = 0; j < K; j++) G[j][i] = RRE_H[i - K][j];
        }
    }

    // swap back and get G
    for(int i = 0; i < N; i++){
        int cur = i;
        int next = change[cur];
        while(next != cur){
            for(int j = 0; j < K; j++) swap(G[j][cur], G[j][next]);
            change[cur] = change[next];
            change[next] = next;
            next = change[cur];
        }
    }

    // HG^T = 0
    int error = 0;
    for(int i = 0; i < M; i++){
        for(int j = 0; j < K; j++){
            int tmp = 0;
            for(int k: graph.CHK[i].neighbor) tmp ^= G[j][k];
            if(tmp != 0) {
                cout << "ERROR!!! " << error++ << endl;
                break;
            }
        }
    }

    outFile << N << " " << K << " " << M << " " << dc << " " << dv << endl;
    for(int i = 0; i < K; i++){
        for(int j = 0; j < N; j++){
            outFile << G[i][j] << " ";
        }
        outFile << endl;
    }

    for(int i = 0; i < M; i++){
        for(int j: graph.CHK[i].neighbor){
            outFile << j << " ";
        }
        outFile << endl;
    }

    //outFile.close();

    return 0;
}