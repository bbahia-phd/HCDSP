function W = Kaiser_sinc(t,N,a);
% Calculate weights via Kaiser window tappered sinc
W = besseli(0,a*sqrt(1-(t/N).^2))/besseli(0,a) .* sinc(t);
end