#!/usr/bin/bash

#
# Run dilution series simulations in parallel
#

NSIM=$1

# no error first
for N in 17 35 75; do
	for IU in 0.2 2 10 20 50; do
		Rscript --vanilla dilutionSimRun.R $N $IU $NSIM ne 2>&1 &
		sleep 1
	done
done

wait $( ps | grep R | awk '// {print $1}' )

for N in 17 35 75; do
	for IU in 0.2 2 10 20 50; do
		Rscript --vanilla dilutionSimRun.R $N $IU $NSIM sf 2>&1 &
		sleep 1
	done
done

wait $( ps | grep R | awk '// {print $1}' )

for N in 17 35 75; do
	for IU in 0.2 2 10 20 50; do
		Rscript --vanilla dilutionSimRun.R $N $IU $NSIM df 2>&1 &
		sleep 1
	done
done

wait $( ps | grep R | awk '// {print $1}' )

