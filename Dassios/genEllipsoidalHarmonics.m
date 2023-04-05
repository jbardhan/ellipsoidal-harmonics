function V = genEllipsoidalHarmonics(ellipseGrid)

np = prod(size(ellipseGrid.xs));

maxOrder = 2; % starts from 0
V1 =makeHarmInterpolants(maxOrder,ellipseGrid,ellipseGrid.xs,ellipseGrid.ys,ellipseGrid.zs);
V2 =makeHarmInterpolants(maxOrder,ellipseGrid,ellipseGrid.xs,-ellipseGrid.ys,ellipseGrid.zs);
V3 =makeHarmInterpolants(maxOrder,ellipseGrid,ellipseGrid.xs,-ellipseGrid.ys,-ellipseGrid.zs);
V4 =makeHarmInterpolants(maxOrder,ellipseGrid,ellipseGrid.xs,ellipseGrid.ys,-ellipseGrid.zs);

V =[V1;V2;V3;V4];
return;

x1=p(1);
x2=p(2);
x3=p(3);
h1=h(1);
h2=h(2);
h3=h(3);

aR = -(x1^2+x2^2+x3^2);
bR = x1^2*(h2^2+h3^2)+h2^2*x2^2+h3^2*x3^2;
cR = -1 - x1^2*h2^2*h3^2;
ff=roots([aR bR cR]);

aR2 = -aR;
bR2 = -bR;
cR2 = -1 + x2^2*h2^2*h3^2;
ff2=roots([aR2 bR2 cR2]);
