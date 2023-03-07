function plotPerf_bo(report)
%
%
%

% Number of iterations and elapsed time
figure(14)
%plotWellSols(wellSols, cumsum(schedule.step.val))
%%
timing = getReportTimings(report);
subplot(2,1,1);
plot([timing.Iterations], '.-')
xlabel('Step')
ylabel('Iterations / step')
%%
subplot(2,1,2);
plot(cumsum([timing.Total])/hour)
xlabel('Step')
ylabel('Elapsed runtime [h]')

end