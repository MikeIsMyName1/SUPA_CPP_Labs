// Michael McFadden
// 09/12

#include "FiniteFunctions.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

int main(){

    // Having studied the source code of the FiniteFunction class, I think I now understand it
    // I need to define an instance with appropriate parameters
    // Then I need to call plotData and plotFunction
    // Then I need to destroy the instance to activate generatePlot

    // But first, I need to import the data from MysteryData15102.txt as a vector of values
    std::ifstream dataInput;
    dataInput.open("MysteryData15102.txt");

    // Check to see if file opened successfully
    if (dataInput.fail()){
        std::cout << "Could not open file." << std::endl;
        exit(1);
    }

    // Now we extract the data to a vector
    double value;
    int numLines=10000;
    std::vector<double> randData;
    randData.reserve(numLines);

    for (int i=1; i<=numLines; i++){
        dataInput >> value;
        randData.push_back(value);
    }

    dataInput.close();

    // Now I can define and use the class instance

    FiniteFunction fin(-10,10,"finFile");
    fin.plotData(randData,50,true);
    fin.plotFunction();
    fin.plotSynth();

    cauchyDist cauchy(-10,10,"cauchyFile");
    cauchy.plotData(randData,50,true);
    cauchy.plotFunction();
    cauchy.plotSynth();

    NormDist norm(-10,10,"normFile");
    norm.plotData(randData,50,true);
    norm.plotFunction();
    norm.plotSynth();

    crystalDist crystal(-10,10,"crystalFile");
    crystal.plotData(randData,50,true);
    crystal.plotFunction();
    crystal.plotSynth();
    return 0;
}

// After playing around with parameters (apologies for having them in the function file, they should
// really be in this test file being called by the class instances), I have concluded that the Cauchy-
// Lorentz distribution best fits the data supplied, and so is likely responsible for the distribution.

// I have managed to get the plots working for all functions and supplied data, before implementing a unilateral synthetic data function in the base class - for some reason,
// this works with the original inverse square function and the Cauchy-Lorentz distribution, but not with the negative Crystal Ball distribution or the normal distribution.
// The only thing I can see in common with these two functions is that they both contain exponentials. Strangely, when I tried running Test without plotSynth in the NormDist
// or crystalDist code sections, an error still occurred, suggesting that I have somehow influenced the rest of the classes' functions when adding my new plot, breaking them.
// I cannot figure it out, so I will have to leave it be - maybe you can tell me what has gone wrong if you know. I have attached the old working plots for those two classes.