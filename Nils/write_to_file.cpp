#include <iostream>
#include <string>
#include <fstream>

int main(){

    // set filename
    std::string filename = "output.txt";
    // create and open the output file.
    // More correctly: creating "output file stream"
    // and connecting it to the filename
    std::ofstream ofile;
    // the following opens and overwrites
    ofile.open(filename);

    // to append use:
    // ofile.open(filename, std::ofstream::app);

    // send som text to the file
    ofile << "Some output to out into file." << std::endl;

    // close the output file
    ofile.close();

    return 0;
}