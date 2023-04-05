function integrands = integrandForExterior(t, n, p, a, b, c)

h = sqrt(a^2 - b^2);
k = sqrt(a^2 - c^2);
lameVals = calcLame(1./t, n, p, a, b, c)';
%integrands = 1./(lameVals.^2 .* sqrt(t.^2 - h^2) .* sqrt(t.^2 - k^2));
integrands = 1./(lameVals.^2 .* sqrt(1-k^2*t.^2) .* sqrt(1-h^2*t.^2));
%plot(t,integrands,'.'); hold on;