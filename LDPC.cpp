#include <iostream>
#include <stdlib.h>
#include "RandVar.h"
#include "Graph.h"
#include <fstream>
#include <vector>
#include <string>

// Hi!!!

using namespace std;

vector<vector<int>> H;
vector<vector<int>> RRE_H;
vector<vector<int>> G;
vector<int> information;
vector<int> codeword;
vector<int> result;
vector<long double> code;
vector<long double> pre_sum;
vector<long double> post_sum;

//double SNR[3] = {1.5, 1.75, 2.0};                              // N = 2640 K = 1320 M = 1320 
//double SNR[6] = {0.9844, 1.289, 1.584, 1.868, 2.411, 3.046};   // N = 816 K = 410 M = 406                 
//double SNR[3] = {0.8279, 1.214, 1.584}; // N = 8000 K = 4000 M = 4000
double SNR[10] = {1.0, 1.2, 1.4, 1.5};

double variable_node_operation(double a, double b);
double check_node_operation(double a, double b);
double zech_logarithm(double x);

int main(int argc, char* argv[]){
    int K, N, M, dc1, dc2;
    int position;
    int index;
    int threshold;
    int errblk_1, errblk_2;
    long long int round;
    bool correct = false;
    long double bit_ep1 = 0;
    long double bit_ep2 = 0;
    double R;
    ifstream inFile;
    ofstream outFile;

    // Read information of parity-check matrix
    inFile.open("database_QC_QAM.txt");
    //cout << "SNR = 1.8 Normal" << endl;
    inFile >> N >> K >> M >> dc >> dv;

    index = N - K;
    R = 1.0 * K / N;
    Graph graph(N, M);

    // Read in Parity-Check Matrix
    G.assign(K, vector<int>(N, 0));
    for(int i = 0; i < K; i++){
        for(int j = 0; j < N; j++){
            inFile >> G[i][j];
        }
    }

    H.assign(M, vector<int>(N, 0));
    for(int i = 0; i < M; i++){
        for(int j = 0; j < dc; j++){
            inFile >> position;
            graph.CHK[i].neighbor.push_back(position);
            graph.VAR[position].neighbor.push_back(i);
            H[i][position] = 1;
        }
    }

    /*int error = 0;
    for(int i = 0; i < M; i++){
        for(int j = 0; j < K; j++){
            int tmp = 0;
            for(int k: graph.CHK[i].neighbor) tmp ^= G[j][k];
            if(tmp != 0) {
                cout << "ERROR!!! " << i << endl;
                break;
            }
        }
    }*/

    code.assign(N, 0);
    codeword.assign(N, 0);
    result.assign(N, 0);
    information.assign(K, 0);
    pre_sum.assign(N + 1, 0);
    post_sum.assign(N + 1, 0);

    //cout << "Eb / N0(dB)               STD_DEV                  bit_ep1                  bit_ep2\n";

    for(int testcase = 0; testcase < 4; testcase++){
        cout << "SNR = " << SNR[testcase] << "(dB)" << endl;
        bit_ep1 = 0;
        information[K - 6] = 0;
        for(int i = 1; i <= 5; i++) information[K - i] = 1;

        errblk_1 = 100;
        round = 0;
        while(errblk_1 > 0){
            correct = false;
            for(int i = 0; i < K; i++) information[i] = information[(i - 6 + K) % K] ^ information[(i - 5 + K) % K];
            for(int i = 0; i < N; i++){
                for(int j = 0; j < K; j++) codeword[i] ^= information[j] * G[j][i];
            }
        
    //////////////////////////////////////////////////////////////////////////////  Noise  //////////////////////////////////////////////////////////////////////////////
            double n1 = 0;
            double n2 = 0;
            STD_DEV = sqrt(pow(10, -SNR[testcase] / 10) / (2 * R));

            // adding noise, y = x + z
            for(int i = 0; i < N; i++) {
                if(i % 2 == 0) {
                    Normal(n1, n2, STD_DEV);
                    if(codeword[i] == 0) code[i] = 1;
                    else code[i] = -1;
                    code[i] += n1;
                }
                else {
                    if(codeword[i] == 0) code[i] = 1;
                    else code[i] = -1;
                    code[i] += n2;
                }
            }
            //for(int i = 0; i < N; i++) cout << codeword[i];
            //cout << endl;

    ////////////////////////////////////////////////////////////////////////////// Decoding /////////////////////////////////////////////////////////////////////////////

        ////// initialize ////////
            for(int i = 0; i < N; i++){
                graph.VAR[i].INIT = 2 * code[i] / (STD_DEV * STD_DEV);
                for(int j: graph.VAR[i].neighbor) graph.VAR[i].OUT[j] = graph.VAR[i].INIT;
            }

            threshold = 50;
            while (!correct && threshold > 0){

                correct = true;

                ////// Bottom-Up /////////
                for(int i = 0; i < M; i++){
                    // calculate prefix and postfix sum
                    pre_sum[0] = graph.VAR[graph.CHK[i].neighbor[0]].OUT[i];
                    post_sum[0] = graph.VAR[graph.CHK[i].neighbor[dc - 1]].OUT[i];
                    for(int j = 1; j < dc; j++) pre_sum[j] = check_node_operation(pre_sum[j - 1], graph.VAR[graph.CHK[i].neighbor[j]].OUT[i]);
                    for(int j = 1; j < dc; j++) post_sum[j] = check_node_operation(post_sum[j - 1], graph.VAR[graph.CHK[i].neighbor[dc - j - 1]].OUT[i]);

                    // calculate message
                    graph.CHK[i].OUT[graph.CHK[i].neighbor[0]] = post_sum[dc - 2];
                    graph.CHK[i].OUT[graph.CHK[i].neighbor[dc - 1]] = pre_sum[dc - 2];
                    for(int j = 1; j <= dc - 2; j++) graph.CHK[i].OUT[graph.CHK[i].neighbor[j]] = check_node_operation(pre_sum[j - 1], post_sum[dc - 2 - j]); 
                }

                ////// Top-Down and Decision //////////
                for(int i = 0; i < N; i++){
                    double tmp = graph.VAR[i].INIT;
                    for(int j: graph.VAR[i].neighbor) tmp += graph.CHK[j].OUT[i];
                    for(int j: graph.VAR[i].neighbor) graph.VAR[i].OUT[j] = tmp - graph.CHK[j].OUT[i];
                    graph.VAR[i].LLR = tmp;

                    if(tmp >= 0) result[i] = 0;
                    else result[i] = 1;
                }

                for(int i = 0; i < M && correct; i++){
                    int tmp = 0;
                    for(int j: graph.CHK[i].neighbor) tmp ^= result[j];
                    if(tmp != 0) correct = false;
                }

                threshold--;
            }

            for(int i = 0; i < N / 2; i++){
                if(result[i] != codeword[i]){
                    errblk_1--;
                    for(int j = i; j < N / 2; j++)
                        if(result[j] != codeword[j]) bit_ep1++;
                    break;
                }
            }
            for(int i = N / 2; i < N; i++){
                if(result[i] != codeword[i]){
                    errblk_2--;
                    for(int j = N / 2; j < N; j++)
                        if(result[j] != codeword[j]) bit_ep2++;
                    break;
                }
            }
            round++;

            if(errblk_1 == 0){
                cout << "BEP 1 : " << round << " " << bit_ep1 / (N * round / 2) << endl;
                errblk_1--;
            }
            if(errblk_2 == 0){
                cout << "BEP 2 : " << round << " " << bit_ep2 / (N * round / 2) << endl;
            }
        }
        //cout << round << " " << bit_ep1 << " " << bit_ep2 << endl;
        //cout << SNR[testcase] << "               " << STD_DEV << "                  " << bit_ep1 / (N * round / 2) << "                  " << bit_ep2 / (N * round / 2) << endl;
    }
    return 0;
}

double variable_node_operation(double a, double b){
    return a + b;
}

double check_node_operation(double a, double b){
    double s1, s2;
    double abs_a, abs_b;
    double p1, p2;
    if(a + b < 0) p1 = -(a + b);
    else p1 = a + b;
    if(a - b < 0) p2 = b - a;
    else p2 = a - b;
    if(a < 0){
        s1 = -1;
        abs_a = -a;
    }
    else {
        s1 = 1;
        abs_a = a;
    }
    if(b < 0){
        s2 = -1;
        abs_b = -b;
    }
    else {
        s2 = 1;
        abs_b = b;
    }
    return s1 * s2 * min(abs_a, abs_b) + zech_logarithm(p1) - zech_logarithm(p2);
}

double zech_logarithm(double x){
    if(x >= 0 && x < 0.196) return 0.65;
    else if(x >= 0.196 && x < 0.433) return 0.55;
    else if(x >= 0.433 && x < 0.71) return 0.45;
    else if(x>= 0.71 && x < 1.05) return 0.35;
    else if(x >= 1.05 && x < 1.508) return 0.25;
    else if(x >= 1.508 && x < 2.252) return 0.15;
    else if(x >= 2.252 && x < 4.5) return 0.05;
    else return 0;
}