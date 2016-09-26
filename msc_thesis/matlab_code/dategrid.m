function [thedates] = dategrid(sd, im, DeltaTime)
    if DeltaTime >= (1/12)
        monthsbetween = months(sd, im)
        nPeriods = monthsbetween/(12*DeltaTime)
        thedates = daysadd(sd,360*DeltaTime*(0:nPeriods),1)
    else    
        daysbetween = daysact(sd, im);
        thedates = daysadd(sd,(0:360*DeltaTime:daysbetween),1);
    end
end

