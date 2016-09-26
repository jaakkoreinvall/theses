function dates = tomatrix(numdates,m,n)
%TOMATRIX Summary of this function goes here
%   Detailed explanation goes here
dates=cell(m,n);
for i=1:m
    for j=1:n
        dates{i,j}=datestr(numdates(i,j));
    end 
end

