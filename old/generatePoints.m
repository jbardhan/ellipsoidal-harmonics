function xyz = generatePoints(nc, pqrFile)
pqrData = readpqr(pqrFile);
xyz = zeros(nc,3);
numatoms = length(pqrData.rads);
for i=1:nc
  whichAtom = -1;
  while whichAtom < 0
	 whichAtom = floor(numatoms*rand(1))+1;
	 if pqrData.rads(whichAtom) < 1e-6
		whichAtom = -1;
	 end
  end
  newpt = pqrData.xyz(whichAtom,:) + ...
			 2*pqrData.rads(whichAtom) * [1 0 0];
  while norm(newpt - pqrData.xyz(whichAtom,:)) > pqrData.rads(whichAtom)
	 newpt = 2*pqrData.rads(whichAtom) * (rand(1,3) - 0.5) + ...
				pqrData.xyz(whichAtom,:);
  end
if 0
  fprintf('%f %f %f is within %f of %f %f %f\n',...
			newpt(1),newpt(2),newpt(3),pqrData.rads(whichAtom),...
			pqrData.xyz(whichAtom,1), ...
			pqrData.xyz(whichAtom,2), ...
			pqrData.xyz(whichAtom,3));
end
  xyz(i,:) = newpt;
end
