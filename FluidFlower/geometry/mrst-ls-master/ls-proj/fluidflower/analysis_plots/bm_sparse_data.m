% Analysis of sparse data and comparison with experiment (table)
% This is done only for the best models, i.e. with 5mm and hysteresis, for
% the 3 models. In the case of model 1, it is done 
clear, close all
% Get sparse data
pth = 'C:/Users/lsalo/matlab/sim_data/mrst/fluidflower/results/';
dir_names = {'benchmark/bmesh5_thickVar_removeS1_modelcase1_D1e-09_Im1_SgMinHyst0.02_kdF3mm_pcD0.3/';
             'benchmark/bmesh5_thickVar_removeS1_modelcase2_D1e-09_Im1_SgMinHyst0.02/';
             'benchmark/bmesh5_thickVar_removeS1_modelcase3_D1e-09_Im1_SgMinHyst0.02/';
             'benchmark/bmesh5_thickVar_removeS1_modelcase3_D3e-09_Im1_SgMinHyst0.02/'};
exp_data = readtable([pth 'experiments/sparse_data']);
m1 = readtable([pth dir_names{1} 'sparse_data']);
m2 = readtable([pth dir_names{2} 'sparse_data']);
m3 = readtable([pth dir_names{3} 'sparse_data']);
m33 = readtable([pth dir_names{4} 'sparse_data']);


% Select values and compute errors with experimental data
e_mean = exp_data.mean;    e_mean(10) = mean(e_mean(10:11)); e_mean(11) = [];
e_std = exp_data.stdDev;   e_std(10) = mean(e_std(10:11)); e_std(11) = [];
e_mean([2:9 11]) = e_mean([2:9 11])*1e3; % g
e_std([2:9 11]) = e_std([2:9 11])*1e3;
m1_mean = m1.p50_mean; m1_mean([2:9 11]) = m1_mean([2:9 11])*1e3; % g 
m2_mean = m2.p50_mean; m2_mean([2:9 11]) = m2_mean([2:9 11])*1e3;
m3_mean = m3.p50_mean; m3_mean([2:9 11]) = m3_mean([2:9 11])*1e3;
m33_mean = m33.p50_mean; m33_mean([2:9 11]) = m33_mean([2:9 11])*1e3;
err1 = (abs(e_mean - m1_mean)./e_mean)*100; err1([3 6 7]) = 0;
err2 = (abs(e_mean - m2_mean)./e_mean)*100; err2([3 6 7]) = 0;
err3 = (abs(e_mean - m3_mean)./e_mean)*100; err3([3 6 7]) = 0;
err33 = (abs(e_mean - m33_mean)./e_mean)*100; err33([3 6 7]) = 0;
err1(10) = 0; err2(10) = 0; err3(10) = 0; err33(10) = 0;
if isinf(err3(3))
    err3(3) = nan;
end
if isinf(err33(3))
    err33(3) = nan;
end
err1_mean = mean(err1, 'omitnan');
err2_mean = mean(err2, 'omitnan');
err3_mean = mean(err3, 'omitnan');
err33_mean = mean(err33, 'omitnan');

% Generate output table
hdrs = {'idx', 'e_mean', 'e_std', 'm1', 'err1', 'm2', 'err2', 'm3', ...
        'err3', 'm33' 'err33'};
vartyp  = cellstr(repmat('double', numel(hdrs), 1));
nr = 12;
sparse = table('Size', [nr, numel(hdrs)], 'VariableTypes', vartyp, ...
               'VariableNames', hdrs);
sparse.idx = {2, '3a', '3b', '3c', '3d', ...
              '4a', '4b', '4c', '4d', 5, 6, 'avg_err'}';
sparse.e_mean =  [e_mean; nan];
sparse.e_std = [e_std; nan];
sparse.m1 = [m1_mean; nan];
sparse.m2 = [m2_mean; nan];
sparse.m3 = [m3_mean; nan];
sparse.m33 = [m33_mean; nan];
sparse.err1 = [err1; err1_mean];
sparse.err2 = [err2; err2_mean];
sparse.err3 = [err3; err3_mean];
sparse.err33 = [err33; err33_mean];