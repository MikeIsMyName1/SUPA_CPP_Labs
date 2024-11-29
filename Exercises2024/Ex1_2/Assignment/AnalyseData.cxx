// Michael McFadden, 29/11/2024

#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
# include "CustomFunctions.h"

// We use the standard namespace, which means that we don't need to prefix our library functions with "std::"
using namespace std;

// We have implemented our functions, so now we can simply request some inputs and run the functions
int main() {

// We first request the filename
string filename;
string message4 = "Please enter the filename you'd like to extract data from.";
Print(message4);
cin >> filename;

// Then we request the error filename
string errorfile;
string message5 = "Please enter the filename you'd like to extract errors from.";
Print(message5);
cin >> errorfile;

// We then use the reading function to create the data vector, then again for the error vector
vector<pair<double, double>> data = readData(filename);
vector<pair<double, double>> errors = readData(errorfile);

// Then we use the squares calculation function to create the squares vector
vector<double> squares = calcSquares(data);

// Then we use the X^Y function for each entry of the vector
for (const auto& point : data) {
        double x = point.first;
        double y = point.second;
        
        // Call the new function to calculate x^y with y rounded
        double result = calcXY(x, y);
        
        std::cout << "x^y: " << result << std::endl;
    }

// We then ask if the user would like to print, and if so, we ask how many lines, and then print
char printChoice;
    string message0 = "Would you like to print the data? (y/n)";
    Print(message0);
    std::cin >> printChoice;

    if (printChoice == 'y' || printChoice == 'Y') {
        // Ask how many lines they would like to print
        int n;
        string message1 = "How many lines would you like to print?";
        Print(message1);
        std::cin >> n;

        Print(data, squares, n);

    }

// We use the line fitting function to find m and c for the data
pair<pair<double, double>,pair<double, double>> lineParams = lineFit(data, errors);
    double m = lineParams.first.first;
    double c = lineParams.first.second;
    double chiSquare = lineParams.second.first;
    int ndf = lineParams.second.second;

    // We then print the fitted line
    Print(lineParams);

// We open the file that we are putting the best fit line into
ofstream outFile("best_fit_line.txt");
    if (outFile.is_open()) {
        // We create a string for the equation of the line
        stringstream equation;
        equation << "y = " << m << "x + " << c << ", chi squared = " << chiSquare << ", NDF = " << ndf;

        // We then write the equation string to the file and announce that we've done so
        outFile << equation.str() << endl;
        string message2 = "Best fit equation saved to 'best_fit_line.txt'.";
        Print(message2);

        outFile.close();

        // If we can't open the file, we display an error message
    } else {
        string message3 = "Unable to open file for writing!";
        Print(message3);

    }

return 0; 

}

// It works!