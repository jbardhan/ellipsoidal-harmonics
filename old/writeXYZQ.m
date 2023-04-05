function writeXYZQ(filename, xyz, q)

if size(xyz,1) ~= length(q)
  fprintf('ERROR: xyz and q vectors have different dimensions!\n');
  return;
end

INFILE = fopen(filename, 'w');
for i=1:length(q)
  fprintf(INFILE, '%f  %f  %f  %f\n',xyz(i,1),xyz(i,2),xyz(i,3),q(i));
end
fclose(INFILE);
