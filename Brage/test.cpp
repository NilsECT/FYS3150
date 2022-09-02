
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

int main() {
    std::vector<int> tall(10);

    for (int t : tall) {
        std::cout << t << std::endl;
    }

    return 0;
}