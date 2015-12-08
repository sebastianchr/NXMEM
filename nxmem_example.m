%PbTe %Uiso=0.0053
clear all
close all
input_structure=[];
input_structure.element={'Pb' 'Te'};
input_structure.verbose=3;
input_structure.outlier_limit=5;
input_structure.sigma_reject=3;
input_structure.cell=[6.435489  0 0 ; 0 6.435489  0 ; 0 0 6.435489];

%Use isotropic model to deconvolve
input_structure.filebase='pbte_oc0_powder';
input_structure.outputfile='pbte_oc0_proffit_NXMEM.BayMEM';
NXMEM2('pbte_oc0_proffit.BayMEM', input_structure);
NXMEM2('pbte_oc0_proffit.BayMEM','pbte_oc0_powder',{'Pb' 'Te'});






