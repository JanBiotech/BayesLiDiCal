/*
 * Copyright (c) 2019 Anthony J. Greenberg
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


/// R interface to the MCMC sampler
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 1.0
 *
 * Contains the R interface to quantile limited dilution essay model fitting.
 *
 */

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <Rcpp.h>

#include "model.hpp"

//' Run MCMC
//'
//' Fits a binomial model to well count data. Uses Metropolis-Hastings MCMC to generate samples from the posterior distribution of IUPM.
//'
//' @param nPos number of positive wells
//' @param nWells total number of wells for each dilution
//' @param dilFrac dilution fractions
//' @param nChains number of chains
//' @param nBurnin number of burn-in iterations
//' @param nSample number of sampling iterations
//'
//' @return A list containing chains of IUPM values, a vector of chain IDs, and a vector of acceptance rates
//'
//' @export
//[[Rcpp::export]]
Rcpp::List runSampler(const std::vector<double> &nPos, const std::vector<double> &nWells, const std::vector<double> &dilFrac, const int32_t &nChains, const int32_t &nBurnin, const int32_t &nSample){
	if (nBurnin <= 0) {
		Rcpp::stop("ERROR: number of burn-in steps must be positive");
	}
	if (nSample <= 0) {
		Rcpp::stop("ERROR: number of sampling steps must be positive");
	}
	if (nChains <= 0) {
		Rcpp::stop("ERROR: number of chains must be positive");
	}
	uint32_t Nb = static_cast<uint32_t>(nBurnin);
	uint32_t Ns = static_cast<uint32_t>(nSample);
	try {
		std::vector<double> iupm;
		std::vector<uint32_t> accept;
		std::vector<int32_t> chainID;
		for (int32_t c = 1; c <= nChains; ++c) {
			BayesicSpace::BayesQLD qld(nPos, nWells, dilFrac);
			qld.sampler(Nb, Ns, iupm, accept);
			for (uint32_t i = 0; i < Ns; ++i) {
				chainID.push_back(c);
			}
		}
		return Rcpp::List::create(Rcpp::Named("iupm", iupm), Rcpp::Named("chainID", chainID), Rcpp::Named("acceptance", accept));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}
	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}


