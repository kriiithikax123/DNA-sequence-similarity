#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include <chrono>

#include "jaccard.hpp"

using namespace std;

string generateRandomDNA(int length, mt19937& generator){
    const string bases="ACGT";

    uniform_int_distribution<int> distribution(0,3);

    string dna;
    dna.reserve(length);

    for(int i=0;i<length;i++){
        dna+=bases[distribution(generator)];
    }

    return dna;
}

void mutateDNA(string& dna,mt19937& generator){
    uniform_int_distribution<int> positionDistribution(0,dna.length()-1);

    uniform_int_distribution<int> baseDistribution(0,3);

    const string bases="ACGT";

    int position=positionDistribution(generator);
    char oldBase=dna[position];

    char newBase=oldBase;

    while(newBase==oldBase){
        newBase=bases[baseDistribution(generator)];
    }

    dna[position]=newBase;
}


int main(){
    random_device rd;
    mt19937 generator(rd());

    ofstream file("results/similarity.csv");

    if(!file.is_open()){
        cerr<<"Error: Could not open CSV file."<<endl;
        return 1;
    }

    file<<"sequence_length,mutations,k,similarity,runtime_microseconds\n";

    int lengths[]={100,1000,10000};

    for(int length:lengths){
        string A=generateRandomDNA(length,generator);
        string B=A;

        for(int mutations=0;mutations<=10;mutations++){
            for(int k=1;k<=32;k++){
                auto countsA=countKmers(A,k);
                auto countsB=countKmers(B,k);

                const int repetitions=100;

                double similarity=0.0;
                double totalRuntime=0.0;

                for(int repetition=0;repetition<repetitions;repetition++){
                    auto start=chrono::steady_clock::now();

                    similarity=calculateJaccard(countsA,countsB);

                    auto end=chrono::steady_clock::now();

                    totalRuntime+=chrono::duration<double,micro>(end-start).count();
                }

                double runtime=totalRuntime/repetitions;

                file<<length<<","<<mutations<<","<<k<<","<<similarity<<","<<runtime<<"\n";
            }
            cout<<"Finished mutations: "<<mutations<<endl;

            if(mutations<10){
                mutateDNA(B,generator);
            }
        }

        cout<<"Finished sequence length: "<<length<<endl;
    }

    file.close();

    cout<<"Results saved to results/similarity.csv"<<endl;

    return 0;

}