function A=getcell_BayMEM(inputfile)

%read BayMEM file
fid = fopen(inputfile);
if fid==-1; error('Inputfile not found'); return; end
x=fread(fid,'*char')';
fclose(fid);

aa=regexp(x, '(?-s)(?m)cell[\t ]+([0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+[\t ]+[0-9\.]+)','tokens');
cellstr=aa{1}{1};
cell=str2num(cellstr);

%from giacovazzo s 75
A=cell2A(cell);
%first row in A is the first lattice vector and so forth.
%Be carefull som may use the transpose definition
end