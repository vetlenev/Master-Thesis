#!/bin/bash

#mkdir mrst-dev
#cd mrst-dev
cd /Users/lluis/Documents/MATLAB/mrst-dev/
repos='core autodiff multiscale solvers model-io visualization co2lab'
for r in $repos; do
  git clone https://lsalo@bitbucket.org/mrst/"mrst-${r}.git"
done
cp ~/matlab/mrst-dev/checkout-mrst/startup_user.m mrst-core/

