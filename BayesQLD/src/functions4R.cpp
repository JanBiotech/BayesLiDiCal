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


///
/** \file
 * \author Anthony J. Greenberg
 * \copyright Copyright (c) 2019 Anthony J. Greenberg
 * \version 1.0
 *
 *
 */

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <Rcpp.h>

#include "model.hpp"

//[[Rcpp::export]]
Rcpp::List testLP(const std::vector<double> &nPos, const std::vector<double> &nWells, const std::vector<double> &dilFrac, const int32_t &nBurnin, const int32_t &nSample){
	try {
		BayesicSpace::BayesQLD qld(nPos, nWells, dilFrac);
		std::vector<double> iupm;
		std::vector<uint32_t> accept;
		if (nBurnin <= 0) {
			Rcpp::stop("ERROR: number of burn-in steps must be positive");
		}
		if (nSample <= 0) {
			Rcpp::stop("ERROR: number of sampling steps must be positive");
		}
		uint32_t Nb = static_cast<uint32_t>(nBurnin);
		uint32_t Ns = static_cast<uint32_t>(nSample);
		qld.sampler(Nb, Ns, iupm, accept);
		return Rcpp::List::create(Rcpp::Named("iupm", iupm), Rcpp::Named("acceptance", accept));
	} catch(std::string problem) {
		Rcpp::stop(problem);
	}
	return Rcpp::List::create(Rcpp::Named("error", "NaN"));
}


