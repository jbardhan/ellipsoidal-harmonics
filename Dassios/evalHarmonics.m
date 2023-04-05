function harms = evalHarmonics(maxOrder,ellipse,MU,NU,verts)

mu = mean(MU(verts));
nu = mean(NU(verts));
harmsmu = evalLame(mu,maxOrder,ellipse);
harmsnu = evalLame(nu,maxOrder,ellipse);
harms = harmsmu.*harmsnu;

