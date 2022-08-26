
#include "test.hpp"
#include <string>
#include <fstream>

int main(int argc, char* argv[]) {

    if (argc != 4) {
        std::string filename = argv[0];

        std::cerr << "Expected 3 arguments, got " << argc-1 << std::endl;
        
        return 1;
    }
    // Print a message to the screen
    int num1 = atoi(argv[1]);
    int num2 = atoi(argv[2]);

    std::string filename = argv[3];

    int added = add(num1, num2);

    std::string num1s = std::to_string(num1), num2s = std::to_string(num2), addeds = std::to_string(added);

    std::string answer = num1s + " + " + num2s + " = " + addeds;

    std::ofstream ofile;
    ofile.open(filename);

    ofile << answer << std::endl;

    ofile.close();

    // Nice!
    return 0;
}