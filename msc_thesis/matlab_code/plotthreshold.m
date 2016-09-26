function [] = plotthreshold(cva, dva, cvac, dvac, cvap, dvap, cvacpt0, ...
    dvacpt0, cvacpt1000, dvacpt1000)
figure;
bar([cva dva (cva-dva); cvac dvac (cvac-dvac); cvap dvap (cvap-dvap); cvacpt0 ...
    dvacpt0 (cvacpt0-dvacpt0); cvacpt1000 dvacpt1000 (cvacpt1000-dvacpt1000)]);
Labels = {'No CSAs', 'One-way CSA \newline(collateral posted by CP)', 'One-way CSA \newline(collateral posted by P)', 'Two-way CSA \newline(threshold = 0)', 'Two-way CSA \newline(threshold = 2)'};
set(gca, 'XTick', 1:5, 'XTickLabel', Labels);
legend('CVA','DVA','BCVA','Location','Southeast');
ylabel('Price [%]');
grid on;
end

