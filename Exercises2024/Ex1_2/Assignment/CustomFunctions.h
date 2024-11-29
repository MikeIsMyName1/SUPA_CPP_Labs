// Michael McFadden, 29/11/2024

#ifndef CUSTOMFUNCTIONS_H  // These guard were recommended so I used them
#define CUSTOMFUNCTIONS_H

#include <vector>
#include <string>
#include <utility>

// We declare the functions
std::vector<std::pair<double, double>> readData(const std::string& filename);
std::vector<double> calcSquares(const std::vector<std::pair<double, double>>& data);
std::pair<std::pair<double, double>,std::pair<double, double>> lineFit(const std::vector<std::pair<double, double>>& data, const std::vector<std::pair<double, double>>& errors);

// These are our overloaded print functions; for the data vector, line parameter pair, and message strings, respectively
void Print(const std::vector<std::pair<double, double>>& data, const std::vector<double>& squares, int n);
void Print(const std::pair<std::pair<double, double>,std::pair<double, double>>& lineParams);
void Print(const std::string& message);

double calcXY(double x, double y);

#endif