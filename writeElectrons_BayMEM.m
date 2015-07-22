function writeElectrons_BayMEM(inputfile,outputfile,electrons)

fid = fopen(inputfile);
if fid==-1; disp('Inputfile not found'); return; end
x=fread(fid,'*char')';
fclose(fid);


newx=regexprep(x,'electrons[\t ]*[0-9\.]+',['electrons  ' num2str(electrons)]);

fid = fopen(outputfile,'w');
fprintf(fid,'%s',newx);
fclose(fid);

end
