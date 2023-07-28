function TC = func_TC(T)
t = T/1000;
TC = 1/(7.5408 + 17.962*t + 3.6142*t^2) + 64/t^2.5*exp(-16.35/t); % W/cm/K
end