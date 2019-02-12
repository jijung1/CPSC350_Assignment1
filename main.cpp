/*
Name:                   Jin Jung
Student ID:             2329401
Email:                  jijung@chapman.edu
Course no. & Section:   CS350-02
Assignment:             1
*/

/*
  Driver file to test the functionality of DNAanalysis class.
*/
#include "DNAanalysis.h"

int main(int argc, char** argv)
{
  ofstream ostream;
  ostream.open("jinjung.out");  //clear jinjung.out 
  ostream.close();

  srand(time(NULL));
  string filename = "";
  string filelocation = "";

  if(argc > 1 && argc < 3)
  {
    bool repeat = true;
    filelocation = argv[1];
    if (filelocation.find(".txt") == -1)
    {
      filelocation+=".txt";
    }
    while(repeat)
    {
      ifstream istream(filelocation);
      if(istream.is_open()) //check if file exists
      {
        istream.close();
        DNAanalysis sample1(filelocation);
        sample1.run();
        //once analysis is finish, ask user if another list should be processed.
        string yesno = "";
        cout << "Finished!\nProcess Another list? (enter y/n):";
        cin >> yesno;
        if (cin.fail())
        {
          cout << "Invalid input!\n";
          cin.clear();
          cin.ignore();
          exit(1);
        }
        if ((sample1.toUpperCase(yesno) == "Y") || (sample1.toUpperCase(yesno) == "YES"))
        {
          cout << "Please enter a new file location: ";
          cin >> filelocation;
          if (cin.fail())
          {
            cout << "Invalid input!\n";
            cin.clear();
            cin.ignore();
            exit(1);
          }
          if (filelocation.find(".txt") == -1)
          {
            filelocation+=".txt";
          }

        }
        else
        {
          repeat = false;
        }
      }
      else
      {
        cout << "invalid input! Exiting..\n";
        exit(1);
      }
    }
  }
  else
  {
    cout << "Invalid terminal arguments\nPlease use format: ./[executable] [filelocation]\n";
    return 0;
  }
  getchar();
  return 0;
}
