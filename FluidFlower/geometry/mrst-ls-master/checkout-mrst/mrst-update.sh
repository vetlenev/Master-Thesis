#!/bin/bash

cd ~/matlab/mrst-dev
repos='core autodiff multiscale solvers model-io visualization co2lab'
for r in $repos; do
  cd "mrst-${r}"  
  git pull https://lsalo@bitbucket.org/mrst/"mrst-${r}" master
  cd ..
done
#cd mrst-ls
#git pull mrst-ls master