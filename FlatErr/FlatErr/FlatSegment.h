#pragma once

#include <vector>

typedef std::pair<double, double> dp;

class FlatSegment{
	std::vector<double> pts;
public:
	FlatSegment(double a, double b, size_t n){
		if(!n) return; //temporary solution
		double dx = (b - a)/n;
		pts = std::vector<double>(n+1);
		for(size_t i = 0; i<n+1; ++i) pts[i] = a + dx*i;
	}
	
	FlatSegment(const std::vector<double>& subdivision):pts(subdivision) {}

	dp GetSubsegment(size_t i){
		return std::pair<double, double>(pts[i-1], pts[i]);
	}

	size_t size() {return pts.size() - 1;}
};