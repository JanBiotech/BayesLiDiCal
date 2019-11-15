/*
 * Copyright (c) 2019 JanBiotech, Inc.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/// Model for dilution series
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 JanBiotech, Inc.
 * \version 1.0
 *
 *  Function implementation for estimating the number of positives from a quantal limited dilution assay.
 *
 */

#include <vector>
#include <cmath>

#include "model.hpp"
#include "random.hpp"

using std::vector;
using namespace BayesicSpace;


// BayesQLD methods
BayesQLD::BayesQLD(const vector<double> &pWellN, const vector<double> &totWellN, const vector<double> &dilutionFrac) : posWells_{pWellN}, nWells_{totWellN}, dilution_{dilutionFrac}, lambda_{0.0001} {
	if( !( (pWellN.size() == totWellN.size()) && (totWellN.size() == dilutionFrac.size()) ) ){
		throw string("ERROR: all input vectors must be of the same size");
	}
	theta_ = exp(rng_.rnorm());
};

void BayesQLD::sampler(const uint32_t &Nburnin, const uint32_t &Nsamples, vector<double> &thetaSamp, vector<uint32_t> &accept){
	for (uint32_t iBnin = 1; iBnin <= Nburnin; ++iBnin) { // burn-in phase
		update_();
	}
	for (uint32_t iSamp = 1; iSamp <= Nsamples; ++iSamp) { // sampling phase
		accept.push_back( update_() );
		thetaSamp.push_back(theta_);
	}

}

double BayesQLD::logPost_(const double &theta){
	double lp = 0.0;

	for (size_t i = 0; i < posWells_.size(); ++i) {
		double lnIsum = 0.0;
		for (double di = posWells_[i] + 1.0; di <= nWells_[i]; di += 1.0) {
			lnIsum += log(di);
		}
		double lnKsum = 0.0;
		for (double dk = 2.0; dk <= nWells_[i] - posWells_[i]; dk += 1.0) {
			lnKsum += log(dk);
		}
		lp += lnIsum - lnKsum + posWells_[i]*log( 1.0 - exp(-dilution_[i]*theta) ) + (posWells_[i] - nWells_[i])*dilution_[i]*theta - lambda_*theta;
	}

	return lp;
}

uint32_t BayesQLD::update_(){
	double sd          = 0.4;                                                                  // proposal SD; calibrated to target 64% acceptance
	double lTheta      = log(theta_);
	double lThetaPrime = lTheta + sd*rng_.rnorm();                                             // proposing a move in log-space
	double lAlpha      = logPost_(exp(lThetaPrime)) - logPost_(theta_) + lThetaPrime - lTheta; // MH criterion; the second subtraction accounts for non-symmetry
	double lU          = log(rng_.runifnz());
	if (lU < lAlpha) {
		theta_ = exp(lThetaPrime);
		return 1;
	} else {
		return 0;
	}
}



