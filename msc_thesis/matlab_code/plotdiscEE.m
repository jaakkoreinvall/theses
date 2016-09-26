function [] = plotdiscEE(discEE,SimDates)
% Counterparty discounted EE
figure;
plot(SimDates,discEE)
datetick('x','mmmyy','keeplimits')
title('Discounted Expected Exposure for Each Counterparty');
ylabel('Discounted Exposure ($)')
xlabel('Simulation Dates')
end

