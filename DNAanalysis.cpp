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
#include "DNAanalysis.h"

DNAanalysis::DNAanalysis(string filelocation) //main constructor
{
  this-> filelocation = filelocation;
  this-> numStrands = 0;
  this-> allNucleotides = 0;
  this->allStrands = "";
  ifstream istream(filelocation);
  string temp = "";
  string allStrandsDelim = "";
  bool invalid = false; //to check for invalid strands
  while(!istream.eof()) //while end of file is not reached
  {
    invalid = false;
    getline(istream,temp);
    if (temp.back() == '\r') {  //remove return carriage character from string temp
    temp.pop_back();
    }
    if(temp.length() <= 0)
    {
      cout << "skipping empty line..\n";
      continue;
    }
    temp = this->toUpperCase(temp); //since capitalization is not guaranteed
    for(int i = 0; i < temp.length(); ++i)  //check for valid inputs for all characters within temp
    {
      if(temp[i] != 'A' && temp[i] != 'T' && temp[i] != 'G' && temp[i] != 'C') //only allow legitimate nucleotides
      {
        cout << "invalid strand! Skipping line" << endl;
        invalid = true;
      }
    }
    if (invalid == false) //only include valid inputs to string allNucleotides and string allStrandsDelim
    {
      this->allNucleotides += temp.length();
      this->allStrandsDelim += temp+" ";
      this->allStrands += temp;
      this->numStrands++;
    }
  }
  istream.close(); //all input strands read, close input file stream.
  if (this->numStrands < 2) //minimum of 3 strands required to perform analysis.
  {
    cout << "Your list doesn't have enough valid strands. Exiting..\n";
    exit(1);
  }
}

int DNAanalysis::getNumStrands() //function returns number of strands read from input file
{
  return this->numStrands;
}
int DNAanalysis::lenSum() // function returns sum of all nucleotides from input file
{
    return this-> allNucleotides;
}

double DNAanalysis::lenMean() //function returns lenSum()/numStrands(). need to cast to type double
{
  return (static_cast<double>(this->allNucleotides) / this->numStrands);
}
double DNAanalysis::lenVar() //function returns lenStdDev()*lenStdDev();
{
  return(pow(this->lenStdDev(),2.0f));
}
double DNAanalysis::lenStdDev() //function returns standard deviation of the length of input dna strands.
{
  //need to calculate the ( sum of (strandlen - lenMean()) ) / numStrands-1
  double diffSum = 0.0f;
  string tempallStrandsDelim = allStrandsDelim;
  while(tempallStrandsDelim.find(char(32)) != std::string::npos)
  {
    int i = 0;
    int templen=0;
    while (tempallStrandsDelim[i] != char(32)) //using whitespace as delimiter
    {
      templen++;
      i++;
    }
    double n1 = templen - this->lenMean();
    double n2 = pow(static_cast<double>(n1),2.0f);
    diffSum += n2;
    tempallStrandsDelim = tempallStrandsDelim.substr(templen+1,tempallStrandsDelim.length()); //get length of strand then remove from tempallStrandsDelim
  }
  return sqrt(diffSum /(this->numStrands - 1));
}
double DNAanalysis::relProb(char a) // function returns relative probability of each nucleotide
{
  int count=0;
  for(int i = 0; i < allStrandsDelim.length(); ++i)
  {
    if(allStrandsDelim[i] == a)
    {
      count++;
    }
  }
  return (static_cast<double>(count) / allNucleotides);
}

double DNAanalysis::bigramProb(char a, char b)//function returns probability of each nucleotide bigram across the entire collection
{
  int numBigrams = allNucleotides / 2;  //assume even number of bigrams
  int count=0;
  for (int i = 0; i < allStrands.length(); i = i+2)
  {
    if(allStrands[i] != static_cast<char>(32) && allStrands[i+1] != static_cast<char>(32))
    {
      if(allStrands[i] == a && allStrands[i+1] == b)
      {
        count++;
      }
    }
  }
  return (static_cast<double>(count) / numBigrams);
}

double DNAanalysis::stdGaussian() //function returns standard gaussian with mean 0 and var 1
{
  double a = (RAND_MAX - rand()) / static_cast<double>(RAND_MAX);
  double b = (RAND_MAX - rand()) / static_cast<double>(RAND_MAX);
  double stdgauss = sqrt(-2*log(a)) * cos(2*M_PI*b);
  return stdgauss;
}
double DNAanalysis::randVar() //function returns a random variable with same mean and variance.
{
  double placeHolder = this->lenStdDev() * this->stdGaussian() + this->lenMean();
  if (placeHolder < 0)
  {
    placeHolder *= -1;
  }
  return placeHolder;
}
void DNAanalysis::writeFile() //function opens a file called yourname.out and prints summary statistics.
{
  ofstream ostream;
  ostream.open("jinjung.out",std::fstream::app);
  ostream << "\n\n------------------------------------------------------------------------------------\n------------------------------------------------------------------------------------\n\n";
  ostream << "Student Name: Jin Jung\nStudent ID: 2329401\n\nSummary statistics:\nLen Sum: " << this->lenSum() << endl; //output sum, mean, var, std dev, relative prob, bigram probs
  ostream << "Len mean: " << this->lenMean() << "\nLen variance: "<< this->lenVar() << "\nLen Standard Dev: "<< this->lenStdDev()<<"\n\n";
  ostream << "Relative prob A: " << this->relProb('A') << "\nRelative prob T: "<< this->relProb('T')<<"\nRelative prob G: "<<this->relProb('G')<<"\nRelative prob C: "<<this->relProb('C')<<endl;
  ostream << "Bigram prob AA: "<<this->bigramProb('A','A')<<" Bigram prob AT: "<<this->bigramProb('A','T')<<" Bigram prob AG: "<<this->bigramProb('A','G')<<" Bigram prob AC: "<<this->bigramProb('A','C')<<endl;
  ostream << "Bigram prob TA: "<<this->bigramProb('T','A')<<" Bigram prob TT: "<<this->bigramProb('T','T')<<" Bigram prob AG: "<<this->bigramProb('T','G')<<" Bigram prob TC: "<<this->bigramProb('T','C')<<endl;
  ostream << "Bigram prob GA: "<<this->bigramProb('G','A')<<" Bigram prob GT: "<<this->bigramProb('G','T')<<" Bigram prob GG: "<<this->bigramProb('G','G')<<" Bigram prob GC: "<<this->bigramProb('G','C')<<endl;
  ostream << "Bigram prob CA: "<<this->bigramProb('C','A')<<" Bigram prob CT: "<<this->bigramProb('C','T')<<" Bigram prob CG: "<<this->bigramProb('C','G')<<" Bigram prob CC: "<<this->bigramProb('C','C')<<endl;
  ostream.close();
}

void DNAanalysis::run()  //in order to calculate the relative frequency of the nucleotides, we need to know the total number of nucleotides for the 1000 generated strands.
            //using this, we can iteratively randomly attach nucleotides and increase the counter of each nucleotide up to its max.
{
  this->writeFile();
  int totalNucleotides = 0, numNucA = 0, numNucT = 0, numNucG = 0, numNucC = 0;
  ofstream ostream;
  ostream.open("jinjung.out",std::fstream::app);
  ostream << "\n\n------------------------------------------------------------------------------------\n------------------------------------------------------------------------------------\n\n";
  ostream << "\n1000 DNA Strands with same mean and variance, and relative probability of nucleotides: \n\n";
  int templength = 0;
  int counter= 0;
  while(counter < 1000)
  {
    templength = static_cast<int>(this->randVar() + .5);
    totalNucleotides += templength;
    if(templength > 0)
    {
      string tempOut = "";

      for(int n = 0; n < templength; n++)
      {
        double tempNum = (RAND_MAX - rand()) / static_cast<double>(RAND_MAX);
        string randomNuc = "";
        if (tempNum < .25)
        {
          if((static_cast<double>(numNucA)/totalNucleotides) < this->relProb('A'))
          {
            randomNuc = "A";
            numNucA++;
          }
          else
          {
            n--;
          }
        }
        else if (tempNum < .5)
        {
          if((static_cast<double>(numNucT)/totalNucleotides) < this->relProb('T'))
          {
            randomNuc = "T";
            numNucT++;
          }
          else
          {
            n--;
          }
        }
        else if(tempNum<.75)
        {
          if((static_cast<double>(numNucG)/totalNucleotides) < this->relProb('G'))
          {
            randomNuc = "G";
            numNucG++;
          }
          else
          {
            n--;
          }
        }
        else
        {
          if((static_cast<double>(numNucC)/totalNucleotides) < this->relProb('C'))
          {
            randomNuc = "C";
            numNucC++;
          }
          else
          {
            n--;
          }
        }
        tempOut += randomNuc;
      }
      ostream << "Strand "<<counter+1 <<": "<< tempOut << endl;
      counter++;
    }
  }
  ostream << "relative probA: " << (static_cast<double>(numNucA)/totalNucleotides) << endl;
  ostream << "relative probT: " << (static_cast<double>(numNucT)/totalNucleotides) << endl;
  ostream << "relative probG: " << (static_cast<double>(numNucG)/totalNucleotides) << endl;
  ostream << "relative probC: " << (static_cast<double>(numNucC)/totalNucleotides) << endl;


  ostream << "\n\n------------------------------------------------------------------------------------\n------------------------------------------------------------------------------------\n\n";
  ostream.close();
}

string DNAanalysis::toUpperCase(string someString)
{
  for(int i = 0; i < someString.length(); ++i)
  {
    if (static_cast<int>(someString[i]) > 96 && static_cast<int>(someString[i]) < 123)
    {
      someString[i] = static_cast<char>(static_cast<int>(someString[i])-32);
    }
  }
  return someString;
}
