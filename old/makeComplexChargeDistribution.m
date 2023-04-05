function [q,xyz]=makeComplexChargeDistribution(c_rad, nc, mindistance, ...
															  l_center, l_rad, l_q, l_xyz);

[new_q, new_xyz] = makeSphereChargeDistribution(c_rad, ...
																mindistance);

inComplexVolume = zeros(length(new_q),1);
ligCenter = [0 0 l_center];
for i=1:length(new_q)
  if (norm(new_xyz(i,:)-ligCenter) > (l_rad+mindistance))
	 inComplexVolume(i) = 1;
  end
end

new_xyz = new_xyz(find(inComplexVolume>0),:);
new_q   = new_q(find(inComplexVolume>0));

picks = randsample(length(new_q),nc);
new_q = new_q(picks);
new_xyz = new_xyz(picks,:);

q = [l_q; new_q];
xyz = [l_xyz+ones(length(l_q),1)*ligCenter; new_xyz];