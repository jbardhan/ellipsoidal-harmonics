function [q, xyz] = makeLigandChargeDistribution(radius, nc,...
																 mindistance)

[q,xyz] = makeSphereChargeDistribution(radius, mindistance);
if (length(q) < nc)
  fprintf('Cannot make %d charges on grid with spacing %f\n',nc,mindistance);
  xyz = 0;
  q = 0;
  return;
end

gridOk = zeros(length(q),1);
for i=1:size(xyz,1)
  if norm(xyz(i,:)) < radius - mindistance
	 gridOk(i) = 1;
  end
end
xyz = xyz(find(gridOk>0),:);
q   = q(find(gridOk>0));

picks = randsample(size(xyz,1),nc);
xyz = xyz(picks,:);
q   = q(picks);