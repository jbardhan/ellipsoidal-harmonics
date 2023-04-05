function [q, xyz] = makeSphereChargeDistribution(radius, ...
																 mindistance)

maxChargeValue = 0.85;

mygrid = -radius:mindistance:radius;
[x,y,z]=meshgrid(mygrid,mygrid,mygrid);

xyz = [reshape(x,numel(x),1) reshape(y,numel(y),1) reshape(z, ...
																  numel(z),1)];
gridOk = zeros(numel(x),1);
for i=1:size(xyz,1)
  if norm(xyz(i,:)) < radius - mindistance
	 gridOk(i) = 1;
  end
end

xyz = xyz(find(gridOk>0),:);
q   = 2*maxChargeValue*(rand(size(xyz,1),1) - 0.5);
