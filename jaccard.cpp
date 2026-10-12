#include "jaccard.hpp"

#include <algorithm>
#include <climits>

using namespace std;

uint64_t encodeBase(char base){
    if(base=='A') return 0;
    if(base=='C') return 1;
    if(base=='G') return 2;
    if(base=='T') return 3;

    return 0;
}

uint64_t encodeKmer(const string& kmer){
    uint64_t value=0;

    for(char base: kmer){
        value=(value<<2) | encodeBase(base);
    }

    return value;
}

uint64_t updateKmer(uint64_t current, char newBase, int k){
    uint64_t mask;

    if(k==32){
        mask=UINT64_MAX;
    }else{
        mask=(1ULL<<(2*k))-1;
    }

    current=(current<<2) & mask;

    current|=encodeBase(newBase);

    return current;
}

unordered_map<uint64_t, uint32_t> countKmers(const string& dna, int k){
    if(k<=0 || k>dna.length()){
        return {};
    }

    unordered_map<uint64_t, uint32_t> counts;

    uint64_t current=encodeKmer(dna.substr(0,k));

    counts[current]++;

    for(int i=k;i<dna.length();i++){
        current=updateKmer(current,dna[i],k);

        counts[current]++;
    }

    return counts;
}

double calculateJaccard(const unordered_map<uint64_t, uint32_t>& countsA, const unordered_map<uint64_t, uint32_t>& countsB){
    uint32_t numerator=0;
    uint32_t denominator=0;

    for(const auto& item: countsA){
        uint64_t kmer=item.first;
        uint32_t countA=item.second;

        uint32_t countB=0;

        if(countsB.find(kmer)!=countsB.end()){
            countB=countsB.at(kmer);
        }

        numerator+=min(countA,countB);
        denominator+=max(countA,countB);
    }

    for(const auto& item: countsB){
        uint64_t kmer=item.first;
        
        if(countsA.find(kmer)==countsA.end()){
            denominator+=item.second;
        }
    }

    if(denominator==0){
        return 1.0;
    }

    return static_cast<double>(numerator)/denominator;
}