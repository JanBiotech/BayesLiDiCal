#
# Copyright (c) 2019 Anthony J. Greenberg
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.
#

#' Estimate IUPM from data
#'
#' Fits a Bayesian binomial model to data, estimating IUPM. The data should be in a data frame that contains at least two columns: (1) positive or negative status and (2) the dilution factor for each well. Positive wells must be 1, negative 0. Each row of the data frame is a well. The dilution factor is the fraction of the initial concentration. If the initial number of cells is 1 million, the result is IUPM. If the initial number is different, however, the reported values must be corrected by dividing the quantiles of the reported value by the initial number of cells.
#'
#' @param data data frame with data
#' @param well.state.colID name or number of the column containing well state data
#' @param dilution.colID name or number of the column with dilution factors
#' @param nBurnin number of burnin iterations (default is 1000, the recommended value)
#' @param nSampling number of sampling iterations per chain (default is 12500, the recommended value)
#' @param nThin chain thinning value (default is 5, the recommended value)
#' @param nChains number of chains (default is 4, the recommended value)
#' @param interval the quantiles of the posterior to report (default is c(0.025, 0.25, 0.5, 0.75, 0.975))
#' @return a list containing the highest-posterior density quantiles of IUPM estimates and MCMC convergence diagnositcs
#'
#' @export
estimateIUPM <- function(data, well.state.colID, dilution.colID, nBurnin=1000, nSampling=12500, nThin=5, nChains=4, interval=c(0.025, 0.25, 0.5, 0.75, 0.975)){
	if (is.na(sum(data[, well.state.colID]))) {
		stop("No missing well state data allowed")
	}
	if (is.na(sum(data[, dilution.colID]))) {
		stop("No missing dilution data allowed")
	}
	dil    <- unique(data[, dilution.colID])
	dilFac <- as.factor(data[, dilution.colID], levels=as.character(dil))
	pos    <- tapply(data[, well.state.colID], dilFac, sum)
	tot    <- tapply(data[, well.state.colID], dilFac, length)
	res    <- runSampler(as.double(pos), as.double(tot), as.double(dil), as.integer(nChains), as.integer(nBurnin), as.integer(nSampling))

	accFrac <- sum(res$acceptance)/length(res$acceptance)
}

