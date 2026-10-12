#ifndef JACCARD_HPP
#define JACCARD_HPP

#include <string>
#include <unordered_map>
#include <cstdint>

using namespace std;

uint64_t encodeBase(char base);

uint64_t encodeKmer(const string& kmer);

uint64_t updateKmer(uint64_t current, char newBase, int k);

unordered_map<uint64_t, uint32_t> countKmers(const string& dna, int k);

double calculateJaccard(const unordered_map<uint64_t, uint32_t>& countsA, const unordered_map<uint64_t, uint32_t>& countsB);    

#endif
