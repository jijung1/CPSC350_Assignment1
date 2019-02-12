/*
Name:                   Jin Jung
Student ID:             2329401
Email:                  jijung@chapman.edu
Course no. & Section:   CS350-02
Assignment:             1
*/

/*
  Class invariant: All objects will have membetr vars int numStrands, int allNucleotides, string allStrandsDelim, string allStrands, and string filelocation.
  Program will read a .txt file specified as an argument in the terminal that contains dna strands and perform statistical analysis on them. 
*/

#ifndef DNAANALYSIS_H
#define DNAANALYSIS_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib> //for rand
#include <ctime> //to seed srand()
#include <string>

using namespace std;

class DNAanalysis
{
  public:
    DNAanalysis(string filelocation); //constructor will open file stream and load all dna strand data to a str mem var
    string toUpperCase(string someString); //prior to concatenating each strand, convert every char to uppercase
    int getNumStrands(); //function returns number of strands read from input file
    int lenSum(); // function returns sum of all nucleotides from input file
    double lenMean(); //function returns lenSum()/numStrands(). need to cast to type double
    double lenVar(); //function returns lenStdDev()*lenStdDev();
    double lenStdDev(); //function returns standard deviation of the length of input dna strands.
    double relProb(char a); // function returns relative probability of each nucleotide
    double bigramProb(char a, char b); //function returns prob ability of each nucleotide bigram across the entire collection
    double stdGaussian(); //function returns standard gaussian with mean 0 and var 1
    double randVar(); //function returns a random variable with same mean and variance.
    void writeFile(); //function opens a file called yourname.out.
                        //writeFile() will then write name, student id, summary statistics.
    void run();
  private:
    int numStrands;
    int allNucleotides;
    string allStrandsDelim;
    string allStrands;
    string filelocation;
};

#endif /*DNAANALYSIS_H*/
