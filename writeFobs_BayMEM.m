function writeFobs_BayMEM(inputfile,outputfile,hkl,A,B,sigma)

fid = fopen(inputfile);
if fid==-1; disp('Inputfile not found'); return; end
x=fread(fid,'*char')';
fclose(fid);
fbegin=regexp(x, '(?-s)(?m)fbegin', 'start');
endf=regexp(x, '(?-s)(?m)endf', 'end');

Fnew= [ sprintf('%s\n','fbegin') sprintf('%4d%4d%4d%14.7f%14.7f%14.7f\n',[hkl A B sigma]') 'endf'];

newx=[x(1:fbegin-1) Fnew x(endf+1:end)];

fid = fopen(outputfile,'w');
fprintf(fid,'%s',newx);
fclose(fid);

end
