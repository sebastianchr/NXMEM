function [output]=adfactor(atom_type, lambda)
%return anomalous scattering factors. Scattering factors are from the
%database used in jana

lambda=lambda(:);

c=2.99792458e8;  %[m/s]
h=4.135669e-18; %[keV s]
%imput lambda in Å
energy_i=h*c./(lambda*1e-10);

atom_type=regexprep(atom_type,'(\<[a-z])','${upper($1)}'); %capitalize first letter
load('jana_ad_database.mat','database')
energy=database.(atom_type).energy; %[keV]
f1=database.(atom_type).f1;
f2=database.(atom_type).f2;

max_energy=max(energy);
if sum(energy_i>max_energy)>0
   error(['X-ray energy exceeds upper limit of ' num2str(max_energy) ' keV']) 
end

%interpolate to get value at energy_i
f1i=interp1(energy,f1,energy_i);
f2i=interp1(energy,f2,energy_i);

output=[lambda energy_i f1i f2i];

end