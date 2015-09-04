#pragma once

#include "FlatSegment.h"
#include "PotentialTrace.h"
#include "EnergySpace.h"

class DiscreteEquation{//smart pointers to be used
	PotentialTrace& F;
	EnergySpace& E;
	FlatSegment& Sg;

	bool mc;
	bool fc;
	std::vector<std::vector<double>> matrix;

public:
	DiscreteEquation(EnergySpace& space, PotentialTrace& trace, FlatSegment& segment):F(trace),E(space),Sg(segment), mc(false), fc(false) {};
	
	void CreateMatrix(){
		std::vector<double> buf = std::vector<double>(Sg.size());
		matrix = std::vector<std::vector<double>>(Sg.size(), buf);

		for(size_t i = 0; i < Sg.size(); ++i){
			DP dp1 = Sg.GetSubsegment(i);
			matrix[i][i] = E(dp1.first, dp1.second, dp1.first, dp1.second);
			for(size_t j = i+1; j < Sg.size(); ++j){
				DP dp2 = Sg.GetSubsegment(j);
				matrix[i][j] = E(dp1.first, dp1.second, dp2.first, dp2.second);
				matrix[j][i] = matrix[i][j];
			}
		}
	}
};