// Michael McFadden
// 09/12

#include <string>
#include <vector>
#include "gnuplot-iostream.h"

#pragma once //Replacement for IFNDEF

class FiniteFunction{

public:
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor
  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  double integral(int Ndiv = 1000); 
  std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  void setOutfile(std::string outfile);
  void plotFunction(); //Plot the function using scanFunction
  
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)

  // Adding two new functions for the Metropolis algorithm into public section of the base class, to be inherited by the custom classes
  void synthesise(); // Will create the synthetic data
  void plotSynth(); // Will automatically call synthesise and then run plotData on the newly updated m_synthData

//Protected members can be accessed by child classes but not users
protected:
  double m_RMin;
  double m_RMax;
  double m_Integral;
  int m_IntDiv = 0; //Number of division for performing integral
  std::string m_FunctionName;
  std::string m_OutData; //Output filename for data
  std::string m_OutPng; //Output filename for plot
  std::vector< std::pair<double,double> > m_data; //input data points to plot
  std::vector< std::pair<double,double> > m_samples; //Holder for randomly sampled data 
  std::vector< std::pair<double,double> > m_function_scan; //holder for data from scanFunction (slight hack needed to plot function in gnuplot)
  bool m_plotfunction = false; //Flag to determine whether to plot function
  bool m_plotdatapoints = false; //Flag to determine whether to plot input data
  bool m_plotsamplepoints = false; //Flag to determine whether to plot sampled data 
  double integrate(int Ndiv);
  std::vector< std::pair<double, double> > makeHist(std::vector<double> &points, int Nbins); //Helper function to turn data points into histogram with Nbins
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 

  // We add a new protected member variable called which will be updated by our "synthesise" function
  std::vector<double> m_synthData;
  
private:
  double invxsquared(double x); //The default functional form
};

// We define an inherited class for the normal distribution

class NormDist : public FiniteFunction{

public:
  NormDist();
  NormDist(double range_min, double range_max, std::string outfile); // Constructors
  double callFunction(double x) override; // Overriding callFunction

private:
  double NormalFunction(double x); // Defining new function for callFunction to call

};

// And similarly for the two other distribution types

class cauchyDist : public FiniteFunction{

public:
  cauchyDist();
  cauchyDist(double range_min, double range_max, std::string outfile); // Constructors
  double callFunction(double x) override; // Overriding callFunction

private:
  double CauchyFunction(double x); // Defining new function for callFunction to call

};


class crystalDist : public FiniteFunction{

public:
  crystalDist();
  crystalDist(double range_min, double range_max, std::string outfile); // Constructors
  double callFunction(double x) override; // Overriding callFunction

private:
  double CrystalFunction(double x); // Defining new function for callFunction to call

};
