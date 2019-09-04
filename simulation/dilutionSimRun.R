#
# Copyright (c) 2019 JanBiotech, Inc.
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

#
# This script runs simulations of dilution series experiments
#
library(compiler)
library(data.table)

args    <- commandArgs(TRUE)
# Number of cells (in millions)
Nml     <- as.integer(args[1])
N       <- Nml*1e6
# IUPM
iupm    <- as.numeric(args[2])
# Number of simulations
Nsim    <- as.integer(args[3])
# simulation type: "ne" for no error, "sf" for sampling first, and "df" for dilution first
simType <- args[4]
n       <- N*iupm/1e6

neSim <- cmpfun(function(N, n){
	pool   <- sample(rep(0:1, times=c(N-n, n)))
	wellID <- factor(rep(c(paste("A", 1:12, sep=""), "B1", "B2", "C1", "C2"),
							times=c(rep(1e6, 12), 2e5, 2e5, 4e4, 4e4)),
					levels=c(paste("A", 1:12, sep=""), "B1", "B2", "C1", "C2"))
	pool   <- pool[1:length(wellID)]
	poolS  <- tapply(pool, wellID, sum)
	res    <- ifelse(poolS == 0, 0, 1)
	return(res)
})

sfSim <- cmpfun(function(N, n){
	pool   <- sample(rep(0:1, times=c(N-n, n)))
	wellID <- NULL
	curN   <- N
	curP   <- 1e6/N
	for (i in 1:12){
		n      <- rbinom(1, curN, curP)
		wellID <- c(wellID, rep(paste("A", i, sep=""), times=n))
		curN   <- curN - n
		curP   <- 1e6/(N - i*1e6) # the proportion is deterministic
	}
	curP   <- 0.550/(N/1e6 - 12)
	bPoolN <- rbinom(1, curN, curP)
	curP   <- 1/2.75
	n      <- rbinom(1, bPoolN, curP)
	wellID <- c(wellID, rep("B1", times=n))
	curP   <- 1/1.75
	bPoolN <- bPoolN - n
	n      <- rbinom(1, bPoolN, curP)
	wellID <- c(wellID, rep("B2", times=n))
	curP   <- 550/750
	bPoolN <- bPoolN - n
	cPoolN <- rbinom(1, bPoolN, curP)
	curP   <- 1/2.75
	n      <- rbinom(1, cPoolN, curP)
	wellID <- c(wellID, rep("C1", times=n))
	cPoolN <- cPoolN - n
	curP   <- 1/1.75
	n      <- rbinom(1, cPoolN, curP)
	wellID <- c(wellID, rep("C2", times=n))
	if(length(pool) > length(wellID)){
		pool <- pool[1:length(wellID)]
	} else if (length(pool) < length(wellID)){
		wellID <- wellID[1:length(pool)]
	}
	wellID <- factor(wellID, levels=c(paste("A", 1:12, sep=""), "B1", "B2", "C1", "C2"))
	poolS  <- tapply(pool, wellID, sum)
	poolS[is.na(poolS)] <- 0
	res    <- ifelse(poolS == 0, 0, 1)
	return(res)
})

dfSim <- cmpfun(function(N, n){
	pool   <- sample(rep(0:1, times=c(N-n, n)))
	wellID <- NULL
	curN   <- N
	curP   <- 0.55e6/N
	bPoolN <- rbinom(1, curN, curP)
	curN   <- curN - bPoolN
	curP   <- 0.55/2.75
	cPoolN <- rbinom(1, bPoolN, curP)
	bPoolN <- bPoolN - cPoolN
	curP   <- 1/2.75
	n      <- rbinom(1, cPoolN, curP)
	wellID <- c(wellID, rep("C1", times=n))
	cPoolN <- cPoolN - n
	curP   <- 1/1.75
	n      <- rbinom(1, cPoolN, curP)
	wellID <- c(wellID, rep("C2", times=n))
	curP   <- 1/2.2
	n      <- rbinom(1, bPoolN, curP)
	wellID <- c(wellID, rep("B1", times=n))
	curP   <- 1/1.2
	bPoolN <- bPoolN - n
	n      <- rbinom(1, bPoolN, curP)
	wellID <- c(wellID, rep("B2", times=n))
	curP   <- 1e6/(N - 0.55e6)
	for (i in 1:12){
		n      <- rbinom(1, curN, curP)
		wellID <- c(wellID, rep(paste("A", i, sep=""), times=n))
		curN   <- curN - n
		curP   <- 1e6/(N - i*1e6) # the proportion is deterministic
	}
	if(length(pool) > length(wellID)){
		pool <- pool[1:length(wellID)]
	} else if (length(pool) < length(wellID)){
		wellID <- wellID[1:length(pool)]
	}
	wellID <- factor(wellID, levels=c(paste("A", 1:12, sep=""), "B1", "B2", "C1", "C2"))
	poolS  <- tapply(pool, wellID, sum)
	poolS[is.na(poolS)] <- 0
	res    <- ifelse(poolS == 0, 0, 1)
	return(res)
})

if (simType == "ne"){
	neSimTable <- data.table(hits=array(replicate(Nsim, neSim(N, n))),
						wellType=rep(c(rep("A", 12), rep("B", 2), rep("C", 2))),
						simulation=rep(1:Nsim, each=16))
	fwrite(neSimTable, file=paste("neSim", Nml, "_", sub("\\.", "", as.character(iupm)), ".tsv", sep=""), quote=FALSE, sep="\t")
} else if (simType == "sf"){
	sfSimTable <- data.table(hits=array(replicate(Nsim, sfSim(N, n))),
						wellType=rep(c(rep("A", 12), rep("B", 2), rep("C", 2))),
						simulation=rep(1:Nsim, each=16))
	fwrite(sfSimTable, file=paste("sfSim", Nml, "_", sub("\\.", "", as.character(iupm)), ".tsv", sep=""), quote=FALSE, sep="\t")
} else if (simType == "df"){
	dfSimTable <- data.table(hits=array(replicate(Nsim, dfSim(N, n))),
						wellType=rep(c(rep("A", 12), rep("B", 2), rep("C", 2))),
						simulation=rep(1:Nsim, each=16))
	fwrite(dfSimTable, file=paste("dfSim", Nml, "_", sub("\\.", "", as.character(iupm)), ".tsv", sep=""), quote=FALSE, sep="\t")
}

