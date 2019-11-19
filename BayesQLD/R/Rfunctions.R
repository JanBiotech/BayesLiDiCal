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
#' @param nBurnin number of burn-in iterations (default is 1000, the recommended value)
#' @param nSampling number of sampling iterations per chain (default is 12500, the recommended value)
#' @param nThin chain thinning value (default is 5, the recommended value)
#' @param nChains number of chains (default is 4, the recommended value)
#' @param probability the probability of the outer HPD quantile to report (default is 0.95)
#' @return a list containing the highest-posterior density quantiles (the outer, inter-quartile range, and mode) of IUPM estimates and MCMC convergence diagnostics
#'
#' @references
#' \insertRef{fazekas82a}{BayesQLD}
#' @export
estimateIUPM <- function(data, well.state.colID, dilution.colID, nBurnin=1000, nSampling=5000, nThin=2, nChains=4, probability=0.95){
	if (is.na(sum(data[, well.state.colID]))) {
		stop("No missing well state data allowed")
	}
	if (is.na(sum(data[, dilution.colID]))) {
		stop("No missing dilution data allowed")
	}
	dil     <- unique(data[, dilution.colID])
	dilFac  <- factor(data[, dilution.colID], levels=as.character(dil))
	pos     <- tapply(data[, well.state.colID], dilFac, sum)
	tot     <- tapply(data[, well.state.colID], dilFac, length)
	res     <- runSampler(as.double(pos), as.double(tot), as.double(dil), as.integer(nChains), as.integer(nBurnin), as.integer(nSampling))
	thinVec <- ifelse(( (1:nSampling)%%nThin ) == 0, TRUE, FALSE)
	chnThn  <- factor(res$chainID[rep(thinVec, nChains)])
	iupmThn <- res$iupm[rep(thinVec, nChains)]
	return(list(acceptance.rate=sum(res$acceptance)/length(res$acceptance), R.hat=split.Rhat(iupmThn, chnThn),
				iupm.HPD=HPDquantile(iupmThn, outr=probability)))
}

#' Find posterior mode
#'
#' Finds the posterior mode from an MCMC sample. If there are more than one, picks one randomly with a warning.
#'
#' @param theta vector of MCMC samples of a variable
#' @return value at the mode
#'
#' @export
pmode <- function(theta){
	dst <- density(theta, adjust = 2)
	mxi <- which(dst$y == max(dst$y))
	if(length(mxi) > 1){
		warning("More than one mode in call to pmode(); picking randomly")
		mxi <- sample(mxi, 1)
	}
	dst$x[mxi]
}

#' HPD interval
#'
#' Finds the highest-posterior density (HPD) interval for the given probability.
#'
#' @param theta vector of MCMC samples
#' @param prob probability (defaults to 0.95)
#' @return parameter values for the lower and upper bound of the HPD
#'
#' @export
HPDint <- function(theta, prob = 0.95){
	nsamp <- length(theta)
	if (nsamp <= 2) stop("vector must have length > 2")

	vals <- sort(theta)
	gap  <- max(1, min(nsamp - 1, round(nsamp * prob)))
	init <- 1:(nsamp - gap)
	mInd <- which.min(vals[init + gap] - vals[init])
	res  <- c(vals[mInd], vals[mInd + gap])
	names(res) <- c("lower", "upper")
	return(res)
}

#' Quantiles of the HPD
#'
#' Finds the mode, the inter-quartile range, and the outer range of the posterior distribution from an MCMC sample. The outer range must be > 0.5.
#'
#' @param theta vector of MCMC samples
#' @param outr outer range (default is 0.95)
#' @return the vector that has the outer range as the first and last elements, the inter-quartile range as the second and penultimate elements, and the mode as the middle element.
#'
#' @export
HPDquantile <- function(theta, outr = 0.95){
	if(outr <= 0.5) stop("Outer margin has to be >= 50%")
	md    <- pmode(theta)
	hpd95 <- HPDint(theta, outr)
	hpd50 <- HPDint(theta, 0.5)
	res   <- c(hpd95[1], hpd50[1], md, hpd50[2], hpd95[2])

	outNm <- paste(c("lower", "upper"), outr*100, sep = "")
	names(res) <- c(outNm[1], "lower50", "mode", "upper50", outNm[2])
	return(res)
}

#' Calculate rank-normalized \eqn{\widehat{R}}
#'
#' Calculates the rank-normalized split-\eqn{\widehat{R}} statistic on regular and folded MCMC samples. Reports the maximum of the two. Folded \eqn{\widehat{R}} captures problems in exploring distribution tails. Values of \eqn{\widehat{R} < 0.01} are considered acceptable.
#'
#' @param theta vector of MCMC samples
#' @param chnFac factor marking chains
#' @return maximum of the folded and unfolded \eqn{\widehat{R}}
#'
#' @references
#' \insertRef{vehtari19a}{BayesQLD}
#' @export
split.Rhat <- function(theta, chnFac){
	S       <- length(theta)
	rnk     <- rank(theta)
	zetaRnk <- rank(abs(theta - median(theta)))
	rn.z    <- qnorm((rnk - 0.5)/S)
	rn.zeta <- qnorm((zetaRnk - 0.5)/S)

	nlev    <- table(chnFac)
	splitLv <- array(rbind(floor(nlev/2), ceiling(nlev/2)))
	splitFc <- factor(paste(chnFac, rep(rep(1:2, nlevels(chnFac)), splitLv), sep=""))
	M       <- nlevels(splitFc)
	N       <- table(splitFc)

	rn.z..    <- mean(rn.z)
	rn.zeta.. <- mean(rn.zeta)
	rn.z.m    <- tapply(rn.z, splitFc, mean)
	rn.zeta.m <- tapply(rn.zeta, splitFc, mean)

	B      <- sum((rn.z.m - rn.z..)^2*N)/(M-1)
	B.zeta <- sum((rn.zeta.m - rn.zeta..)^2*N)/(M-1)
	W      <- sum(tapply((rn.z - rn.z.m[splitFc])^2, splitFc, sum)/(N-1))/M
	W.zeta <- sum(tapply((rn.zeta - rn.zeta.m[splitFc])^2, splitFc, sum)/(N-1))/M

	Nfr       <- mean((N-1)/N)
	N         <- mean(N)
	Rhat.z    <- sqrt(Nfr + B/(W*N))
	Rhat.zeta <- sqrt(Nfr + B.zeta/(W.zeta*N))
	return(max(c(Rhat.z, Rhat.zeta)))
}

