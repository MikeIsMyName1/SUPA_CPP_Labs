// Michael McFadden, 29/11/2024

#include "CustomFunctions.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

// We use the standard namespace, which means that we don't need to prefix our library functions with "std::"
using namespace std;

// We define a function (outside main) for reading the file into a pair-wise vector
vector<pair<double, double>> readData(const string& filename){
    // We create an ifstream object called input2D
    ifstream input2D;

    // We open the file input2D_float.txt into this object
    input2D.open(filename);

    // We will read each line into the variable "line", then process this into x and y via a vector.
    // We initialise the line, vector, and x and y
    string line;
    vector<pair<double, double>> data;
    double x, y;

    // We check that it worked alright
    if(input2D.fail()) {
        cout << "Could not open file." << endl;
        exit(1);
    }

    // We skip over the first line (useless strings) by using the getline function once, outside the loop.
    getline(input2D,line);

    while(getline(input2D, line)){
        // We skip the rest of the loop if the line is blank, as to skip the last line of the file.
        if(line.empty()){
            continue;
        }

        // We now use stringstream to separate the line up into its constituents
        stringstream ss(line);

        // We also discard the comma
        char comma;
        // We check that the stringstream worked and if so, we read the line into a vector pair entry
        if(ss >> x >> comma >> y){
            data.push_back(make_pair(x,y));
        } else {
            cout << "Error processing line:" << line << endl;
        }
    }

    input2D.close();
    return data;
    }

// We define a separate function to calculate the sum of the squares and make a vector of the sums
vector<double> calcSquares(const vector<pair<double, double>>& data){
    vector<double> squares;

    // We take x and y from each entry and add their squares as a new variable, then add that to a vector
    for (const auto& point : data) {
        double x = point.first;
        double y = point.second;
        double sumOfSquares = std::pow(x, 2) + std::pow(y, 2);
        squares.push_back(sumOfSquares);
    }

    return squares;

}

// We define our least-squares fit function
pair<pair<double, double>,pair<double, double>> lineFit(const vector<pair<double, double>>& data, const vector<pair<double, double>>& errors) {
    double sum_x = 0.0, sum_y = 0.0, sum_xx = 0.0, sum_xy = 0.0;
    int N = data.size();

    // We incrementally sum the four different sigma-sum terms that we need for the least-squares fit
    for (const auto& point : data) {
        double x = point.first;
        double y = point.second;
        sum_x += x;
        sum_y += y;
        sum_xx += std::pow(x, 2);
        sum_xy += x * y;
    }

    // Then we use these to calculate the slope and y-intercept using their least-squares formulas
    double m = (N * sum_xy - (sum_x * sum_y)) / (N * sum_xx - pow(sum_x, 2));
    double c = ((sum_xx * sum_y) - (sum_xy * sum_x)) / ((N * sum_xx) - (sum_x * sum_x));

    double chiSquared = 0;
    for (int i = 1; i < N; i++) {
        double x = data[i].first;
        double y = data[i].second;
        double erry = errors[i].second;

        double model = m * x + c;  // The model's predicted value for this x
        chiSquared += pow((y - model), 2) / pow(erry, 2);
    }

    // Degrees of freedom: n - 2 (because we fit two parameters, m and c)
    int ndf = N - 2;

    pair<pair<double, double>,pair<double, double>> fitAndChi = {{m,c},{chiSquared,ndf}};

    // We then return the slope and y-intercept, and chi-squared and NDF, in a pair of pairs
    return fitAndChi;

}


// We define our three different overloads of the print function

// Print function for the x and y data, as well as the sum-of-squares data
void Print(const vector<pair<double, double>>& data, const vector<double>& squares, int n) {
    int totalData = data.size(); // We count the number of rows in our data vector

    if (n > totalData) {
        // We check if n is greater than this number of rows - if it is, we show an error
        cout << "Warning: Requested to print more lines than available. Printing first 5 lines." << endl;
        n = 5;  // Print only the first 5 lines if n is too large
    }

    // We print n lines of the vector, from row 0 to row n, or to the end of the vector if it's less
    for (int i = 0; i < n && i < totalData; ++i) {
        // We extract the first bit of data from each row's pair, then the second, as x and y
        cout << "x: " << data[i].first << ", y: " << data[i].second << ", x^2 + y^2:" << squares[i] << endl;
    }
}

// Print function for the fitted line equation, given the parameters
void Print(const pair<pair<double, double>,pair<double, double>>& lineParams) {
    cout << "Fitted line: y = " << lineParams.first.first << "x + " << lineParams.first.second
    << " . Chi Squared = " << lineParams.second.first << " . NDF = " << lineParams.second.second
    << " ." << endl;

}

// Print function for a string message
void Print(const string& message) {
    cout << message << endl;
}

// A function for finding x^y for each x and y
double calcXY(double x, double y) {
    // We round y to the nearest integer
    y = std::round(y);
    
    double result = pow(x, y);

ofstream outFile("x_to_the_y.txt", std::ios::app);
    if (outFile.is_open()) {
        // We create a string for the equation of the line
        stringstream value;
        value << "X^Y = " << result << endl;

        // We then write the equation string to the file and announce that we've done so
        outFile << value.str() << endl;
        string message7 = "Best fit equation saved to 'x_to_the_y.txt'.";
        Print(message7);

        outFile.close();

        // If we can't open the file, we display an error message
    } else {
        string message3 = "Unable to open file for writing!";
        Print(message3);

    }

    return result; 
}
