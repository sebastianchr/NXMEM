function [hkl A B sigma]=readFobs_BayMEM(inputfile)

fid = fopen(inputfile);
if fid==-1; error('Inputfile not found'); return; end
x=fread(fid,'*char')';
fclose(fid)

fbegin=regexp(x, '(?-s)(?m)fbegin', 'start');
endf=regexp(x, '(?-s)(?m)endf', 'end');

F=x(fbegin:endf);
intsig=regexp(F,['(?-s)([0-9\-]+)[\t ]+([0-9\-]+)[\t ]+([0-9\-]+)(?m)[\t ]+([0-9\.\-]+)[\t ]+([0-9\.\-]+)[\t ]+([0-9\.\-]+)'], 'tokens');

n=length(intsig);
hkl=zeros(n,3); A=zeros(n,1); B=zeros(n,1); sigma=zeros(n,1);

for i =1:n
   temp= str2double(intsig{i});
hkl(i,:)=temp(1:3); A(i)=temp(4); B(i)=temp(5); sigma(i)=temp(6);
end


%format of Fobs in BayMEM file: '%4d%4d%4d%14.7f%14.7f%14.7f'

end