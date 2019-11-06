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
 * Class definition and interface documentation for estimating the number of positives from a quantal limited dilution assay.
 *
 */

#ifndef model_hpp
#define model_hpp

#include <vector>

#include "random.hpp"

using std::vector;

namespace BayesicSpace {
	/** \brief Dilution assay model class
	 *
	 * Keeps the data, fits the model, and saves samples from parameter distributions.
	 */
	class BayesQLD {
	public:
		/** \brief Default constructor */
		BayesQLD() : theta_{0.0} {};
		/** \brief Constructor
		 *
		 * \param[in] pWellN vector with positive well numbers
		 * \param[in] totWellN vector with total well numbers
		 * \param[in] dilutionFrac vector with dilutino fractions
		 * 
		 */
		BayesQLD(const vector<double> &pWellN, const vector<double> &totWellN, const vector<double> &dilutionFrac) : posWells_{pWellN}, nWells_{totWellN}, dilution_{dilutionFrac}, theta_{0.0} {};

		/** \brief Destructor */
		~BayesQLD() {};

		/** \brief Copy constructor
		 *
		 * \param[in] in the object to be copied
		 */
		BayesQLD(const BayesQLD &in) : posWells_{in.posWells_}, nWells_{in.nWells_}, dilution_{in.dilution_}, theta_{in.theta_}, accept_{in.accept_} {};
		/** \brief Copy assignment operator
		 *
		 * \param[in] object to be assigned
		 * \return `BayesQLD` object
		 */
		BayesQLD& operator=(const BayesQLD &in);
		/** \brief Move constructor
		 *
		 * \param[in] in the object to be moved
		 */
		BayesQLD(BayesQLD &&in) : posWells_{move(in.posWells_)}, nWells_{move(in.nWells_)}, dilution_{move(in.dilution_)}, theta_{in.theta_}, accept_{move(in.accept_)} {};
		/** \brief Move assignment operator
		 *
		 * \param[in] object to be assigned
		 * \return `BayesQLD` object
		 */
		BayesQLD& operator=(BayesQLD &&in);

		/** \brief log of the posterior distribution
		 *
		 * \param[in] theta parameter value
		 * \return log-posterior value
		 */
		double logPost(const double &theta) {return logPost_(theta); };
	private:
		/** \brief Number of positive wells at each dilution */
		vector<double> posWells_;
		/** \brief Total number of wells at each dilution */
		vector<double> nWells_;
		/** \brief Dilution fraction */
		vector<double> dilution_;
		/** \brief Current estimate of the number of positives */
		double theta_;
		/** \brief Accept/reject tracking vector */
		vector<uint16_t> accept_;
		/** \brief Sampler from distributions */
		RanDraw rng_;
		/** \brief Log-posterior function
		 *
		 * \param[in] thetaPrime proposed version of \f$\theta\f$
		 * \return Value of the log-posterior
		 */
		double logPost_(const double &thetaPrime);
	};
}


#endif /* model_hpp */
