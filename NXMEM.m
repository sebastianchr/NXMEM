function output=NXMEM(varargin)
%Author: Sebastian Christensen, 25/10/2013 (Grenoble)
%inputfile: BayMEM file
%Note: Not tested for non-centrosymmetric structures
%input_structure:
%Mandatory keywords:
%   element: cell array with name-strings of the elements. Case
%   insensitive.
%   filebase: Partial m80-files (filbase_element.m80)
%Optional keywords:
%   cell: 3x3-matrix defining the unitcell, or 6x1/1x6-matrix = [a b c
%     alpha beta gamma]. If not defined, the cell will be read from the
%     BayMEM file.
%   Outputfile: () file name to write new BayMEM-file to. If not defined, the
%     output will be written to filebase_NXMEM.BayMEM
%   Verbose: (1) Set the verbosity level: 0 = no output from program. 1 =
%     output to command prompt. 3 = text and graphical output.
%   Uiso:  (0.0) scalar. Degree of thermal motion to be deconvolved. Can be
%     negative, to smear the nuclear density.
%   Sigma_reject: (3) remove reflections with intensity significance
%     smaller than this value
%   outlier_limit: (5) if the ratio between F_obs_NDD and F_calc_NDD are
%     larger than this value, the reflection will be removed. It is
%     assumed that the model is good enough that such large values should
%     not occure. Large values are typically encountered for weak
%     reflections.
%   outlier_rule: ('reject') Determines what to do with outliers. 'reject'
%     = outlieing reflections are omitted from calculation. 'calc' = use
%     F_calc_NDD instead of F_obs_NDD
%   df_mean_deviation: (Inf)
%   df_limits: ([-Inf, Inf])
%   weight_exponent: (0)
%
%Largely untested/unfinished keywords:
%   Lebail:
%   ad_correction: performs anomalous dispersion correction based on
%   partial structure factors. To be used in connection with lebail
%     extracted structurefactors.
%
%Latest update: 07/01/2015
%27012015:
%Added the lowerlimit (1/outlier_limit) to the comparison with model
%Added check comparing number of electrons in m80- and
%BayMEM-files. Error can arrise if stochiometry in JANA is incorrect i.e.
%forgot to reset occupancy.

%Dependencies:
%plotplus, readFobs_BayMEM, getcell_BayMEM
%test_variables
output=1;
inputfile=varargin{1};
filebase=varargin{2};
element=varargin{3};
input_structure=struct(varargin{4:end});

%%
%Set standard setting if not defined.
if (isfield(input_structure,'Uiso')); Uiso=input_structure.Uiso;    else Uiso=0; end
if (isfield(input_structure,'sigma_reject'));                       else input_structure.sigma_reject=3; end
if (isfield(input_structure,'weight_exponent'));                    else input_structure.weight_exponent=1; end
if (isfield(input_structure,'outlier_limit'));                      else input_structure.outlier_limit=5; end
if (isfield(input_structure,'df_mean_deviation'));                  else input_structure.df_mean_deviation=Inf; end
if (isfield(input_structure,'df_limits'));                          else input_structure.df_limits=[-Inf Inf]; end
if (isfield(input_structure,'verbose'));                            else input_structure.verbose=3; end
if (isfield(input_structure,'outputfile'));                         else input_structure.outputfile=[filebase '_NXMEM.BayMEM']; end
if (isfield(input_structure,'ad_correction'));                      else input_structure.ad_correction=0; end
if (isfield(input_structure,'outlier_rule'));                       else input_structure.outlier_rule='reject'; end
if (isfield(input_structure,'test'));
    deconvolution_test=input_structure.test;
    if (isfield(input_structure,'deconvolution_factor'));
        deconvolution_factor_type=input_structure.deconvolution_factor;
    else deconvolution_factor_type=3;
    end
else
    input_structure.test=0;
end
%%

if (isfield(input_structure,'cell'));
    %Use the definitions of Giacovazzo (Fundamentals of Crystallography).
    %A is the basis vector matrix. Each row is a separate lattice vector.
    if isequal(size(input_structure.cell),[3,3])
        A=input_structure.cell;
    elseif isequal(size(input_structure.cell),[6,1]) || isequal(size(input_structure.M),[6,1])
        cell=input_structure.cell;
        %from Giacovazzo s 75
        a=cell(1); b=cell(2); c=cell(3);
        alpha=cell(4); beta=cell(5); gamma=cell(6);
        cosalphastar=(cosd(beta)*cosd(gamma)-cosd(alpha))/(sind(beta)*sind(gamma));
        V=a*b*c*(1-cosd(alpha)^2-cosd(beta)^2-cosd(gamma)^2+2*cosd(alpha)*cosd(beta)*cosd(gamma)).^.5;
        cstar=a*b*sind(gamma)/V;
        A=[a 0 0; b*cosd(gamma) b*sind(gamma) 0 ; c*cosd(beta) -c*sind(beta)*cosalphastar 1/cstar];
    else
        error('The format of ''cell'' was not recognized')
    end
    
else
    %get cell from baymem-file
    try
        A=getcell_BayMEM(inputfile);
    catch err
        error('getcell_BayMEM.m for reading cell information from BayMEM file was not found')
    end
    if input_structure.verbose>0
        disp('Using this cell:')
        disp(A)
    end
end

% ###################### End of settings section ###########################
%%
% Initialization. Check that the necessary functions are available.
does_function_exist('plotplus')
does_function_exist('readFobs_BayMEM')
does_function_exist('writeFobs_BayMEM')
does_function_exist('getcell_BayMEM')
does_function_exist('writeElectrons_BayMEM')
does_function_exist('get_f_coefficients')

%%
M=2*pi*(A^-1);

%read Fobs from BayMEM-file:
[hkl_baymem A_obs_EDD B_obs_EDD sigma_obs_EDD]=readFobs_BayMEM(inputfile);

F_obs_EDD=A_obs_EDD+1i*B_obs_EDD;
index_rejection=zeros(size(A_obs_EDD)); %Reflections to be rejected for various reasons

%sigma rejection
index_sigma_reject=reject_based_on_sigma(F_obs_EDD,sigma_obs_EDD,input_structure);
if input_structure.verbose>1
    print_sigma_rejected_reflections(hkl_baymem(index_sigma_reject,:),F_obs_EDD(index_sigma_reject),sigma_obs_EDD(index_sigma_reject))
end
index_rejection=index_rejection | index_sigma_reject;


if (isfield(input_structure,'lebail'));
    if input_structure.lebail
        disp('WARNING: The lebail section has not been tested yet')
        if 0==(isfield(input_structure,'scale_factor'));
            error('Missing scale_factor.');
        end
        
        normF_obs_AD=abs(A_obs_EDD/input_structure.scale_factor);
        sigma_obs_AD=sigma_obs_EDD/input_structure.scale_factor;
        
        %Import partial form factors from m80 files
        [~, A_calc_0_partial,B_calc_0_partial,ff,~, ~]=read_partial_F(filebase,element,M);
        
        %anomalous dispersion correction
        A_calc_0=sum(A_calc_0_partial,2);
        B_calc_0=sum(B_calc_0_partial,2);
        F_calc_0=sqrt(A_calc_0.^2+B_calc_0.^2);
        phase_calc=atan2(B_calc_0,A_calc_0);
        
        %Perform Anomalous dispersion correction - This is not tested. Move
        %to separate function
        if input_structure.ad_correction
            %get ad coefficients from database
            if (isfield(input_structure,'ad'))~=1;
                for i=1:length(element)
                    [~, ~, input_structure.ad.(element{i})(1), input_structure.ad.(element{i})(2)]=adfactor(element{i},input_structure.lambda);
                end
            end
            %AD formfactor
            ff1=zeros(size(ff)); ff2=zeros(size(ff));
            for i=1:length(element)
                ff1(:,i)=ff(:,i)+input_structure.ad.(element{i})(1);
                ff2(:,i)=repmat(input_structure.ad.(element{i})(2),length(ff),1);
            end
            A_calc_AD=sum(A_calc_0_partial./ff.*ff1-B_calc_0_partial./ff.*ff2,2);
            B_calc_AD=sum(B_calc_0_partial./ff.*ff1+A_calc_0_partial./ff.*ff2,2);
            F_calc_AD=sqrt(A_calc_AD.^2+B_calc_AD.^2);
            
            A_obs_AD=cos(phase_calc).*normF_obs_AD;
            B_obs_AD=sin(phase_calc).*normF_obs_AD;
            
            A_obs_0=A_obs_AD-(A_calc_AD-A_calc_0);
            B_obs_0=B_obs_AD-(B_calc_AD-B_calc_0);
            
            %set 000-reflection
            A_obs_0(1)=A_calc_0(1);
            
            F_obs_0=sqrt(A_obs_0.^2+B_obs_0.^2);
            F_calc_0=sqrt(A_calc_0.^2+B_calc_0.^2);
            
            sigma_obs_0=sigma_obs_AD.*F_obs_0./normF_obs_AD;
            
            %Continue with corrected structure factors
            A0=A_obs_0;
            B0=B_obs_0;
            sigma0=sigma_obs_0;
        else % do not do anomalous dispersion correction
            %Phases must be known from model. Intensities from LeBail
            A_obs_AD=cos(phase_calc).*normF_obs_AD;
            B_obs_AD=sin(phase_calc).*normF_obs_AD;
            
            A0=A_obs_AD;
            B0=B_obs_AD;
            sigma0=sigma_obs_AD;
            
            %set 000-reflection
            A0(1)=A_calc_0(1);
        end
        
    end
end
%%
%####################### Deconvolution Section ############################
[F_norm_calc,A_calc_partial,B_calc_partial,ff,Z, hkl_model]=...
    read_partial_F(filebase,element,M);
qnorm_aa=sqrt(sum((M*hkl_model').^2,1))';
m80_electrons=sum(F_norm_calc(1,:));
BayMEM_electrons=F_obs_EDD(1);
xmin=min(qnorm_aa/(4*pi))-0.05*mean(qnorm_aa/(4*pi));
xmax=max(qnorm_aa/(4*pi))+0.05*mean(qnorm_aa/(4*pi));
if abs(m80_electrons-BayMEM_electrons)/(m80_electrons+BayMEM_electrons)>0.001
    error(['Something is wrong with the m80- or BayMEM-inputfiles. The 000-reflection is too different\n m80:' num2str(m80_electrons) '\n BayMEM:' num2str(BayMEM_electrons)])
end




%normalize form factor
if input_structure.verbose>0
    disp('Standard Deconvolution method')
end
ffn=ff./repmat(Z,size(ff,1),1);
A_calc_NDD=sum(A_calc_partial./ffn,2);
A_calc_EDD=sum(A_calc_partial,2);
B_calc_NDD=sum(B_calc_partial./ffn,2);
B_calc_EDD=sum(B_calc_partial,2);

F_calc_NDD=A_calc_NDD+1i*B_calc_NDD;
F_calc_EDD=A_calc_EDD+1i*B_calc_EDD;

is_centrosymmetric=sum(abs(imag(F_calc_EDD)))<1e-10;

%Check ordering of reflections:
[~,ind] = ismember(hkl_baymem,hkl_model,'rows');
F_calc_EDD=F_calc_EDD(ind);
F_calc_NDD=F_calc_NDD(ind);
%ind is the index of the row in hkl_model corresponding to each row in
%hkl_baymem.

df=F_calc_NDD./F_calc_EDD;


if input_structure.test
    switch deconvolution_factor_type
        case 1
            disp('Deconvoluting using method 1')
            ffn=ff./repmat(Z,size(ff,1),1);
            df=1./(sum(F_norm_calc.*ffn,2)./sum(F_norm_calc,2));
        case 2
            disp('Deconvoluting using method 2')
            ffn=ff./repmat(Z,size(ff,1),1);
            df=(sum(F_norm_calc./ffn,2)./sum(F_norm_calc,2));
        case 5
            %deconvolution by Morningstar-Warren approxima tion
            disp('Deconvoluting using method 5 - Morningstar-Warren')
            %stoichiometry averaged
            ffn=ff./repmat(Z,size(ff,1),1);
            weight=A_calc_partial(1,:)./Z;
            weight=weight/sum(weight);
            df=1./(sum(ffn.*repmat(weight,size(ff,1),1),2));
            
            %electron averaged
            %alternative like averaged in patterson paper
            %             weight=A_calc_partial(1,:);
            %             weight=weight/sum(weight);
            %             df=1./(sum(ffn.*repmat(weight,size(ff,1),1),2));
            %
    end
    
end
F_obs_NDD=F_obs_EDD.*df;
sigma_obs_NDD=sigma_obs_EDD.*abs(df);

%remove reflections extremely affected by deconvolution:
index_df_outliers=df < input_structure.df_limits(1) ...
    | df>input_structure.df_limits(2);

index_rejection=index_rejection | index_df_outliers;

%Find estimated nuclear structure factors which deviates
%severely from the value expected for the model structure.
[F_obs_NDD,empiric_weight_correction]=outliers_compared_to_model(F_obs_NDD,F_calc_NDD,input_structure);

index_outliers=isnan(F_obs_NDD);
index_rejection=index_rejection | index_outliers;
additional_outliers=1==(index_outliers-index_sigma_reject);

if input_structure.verbose>0
    fprintf('%.0f additional reflections removed by outlier criterion compared to model\n',sum(additional_outliers))
    F_NDD_deviation=F_obs_NDD./F_calc_NDD;
    print_outlier_reflections(hkl_baymem(index_outliers,:), F_NDD_deviation(index_outliers),input_structure)
end



[~, ind_sort]=sort(qnorm_aa);
df_smooth(ind_sort)=smooth(qnorm_aa(ind_sort),df(ind_sort));
df_smooth(ind_sort)=smooth(qnorm_aa(ind_sort),df_smooth(ind_sort),0.9,'rloess');
df_smooth=df_smooth';
index_outliers_df_mean=abs(df-df_smooth)./df_smooth>input_structure.df_mean_deviation;
additional_outliers=1==(index_outliers_df_mean-index_rejection);
index_rejection=index_rejection | index_outliers_df_mean;

if input_structure.verbose>0
    fprintf('%.0f additional reflections removed as outliers compared to mean deconvolution factor.\n',sum(additional_outliers))
end

ffn=ff./repmat(Z,size(ff,1),1);
weight=A_calc_partial(1,:)./Z;
weight=weight/sum(weight);
df_mw=1./(sum(ffn.*repmat(weight,size(ff,1),1),2));
relative_df=df./df_mw;

if input_structure.verbose > 2
    adata0=[]; adata1=[];
    adata0.hkl=hkl_baymem;
    adata0.sinlambda=qnorm_aa/(4*pi);
    adata1.hkl=hkl_baymem(index_rejection,:);
    adata1.sinlambda=qnorm_aa(index_rejection,:)/(4*pi);
    
    
    min_stol=min(qnorm_aa/(4*pi)); %stol = sinus theta over lambda
    max_stol=max(qnorm_aa/(4*pi));
    if input_structure.df_limits(1)> -Inf
        line([min_stol max_stol], repmat(input_structure.df_limits(1),1,2),'linestyle','--','color','m')
    end
    if input_structure.df_limits(2)< Inf
        line([min_stol max_stol], repmat(input_structure.df_limits(2),1,2),'linestyle','--','color','m')
    end
    
    
    ss = get( 0, 'Screensize' );
    figure;
    if is_centrosymmetric
        set(gcf,'outerposition',[0 ss(4)/2 ss(3)/2 ss(4)/2])
        subplot(1,1,1); hold on;
        h1=plotplus(qnorm_aa/(4*pi), df, adata0,'x');
        
        
        if sum(index_rejection)>0
            plotplus(qnorm_aa(index_rejection)/(4*pi),df(index_rejection),adata1,'ok');
        end
        h2=plot(qnorm_aa(ind_sort)/(4*pi),df_smooth(ind_sort),'-m');
        h3=plot(qnorm_aa(ind_sort)/(4*pi),df_mw(ind_sort),'-g');
        
        
        box on
        ylabel('Deconvolution factor','fontsize',14)
        xlabel('sin(\theta)/\lambda (Å^{-1})','fontsize',14)
        legend([h1 h2 h3], '\it{df}\rm{(\bf{H})}', 'smooth \it{df}', 'MW \it{df}','location', 'northwest')
        xlim([xmin,xmax])
        %             subplot(1,2,2); hold on
        %             plotplus(qnorm_aa/(4*pi),relative_df,adata0,'x');
        %             if sum(index_rejection)>0
        %                 plotplus(qnorm_aa(index_rejection)/(4*pi),relative_df(index_rejection),adata1,'ok');
        %             end
        %             box on
        %             ylabel('df / MW-df','fontsize',14)
        %             xlabel('sin(\theta)/\lambda (Å^{-1})','fontsize',14)
    else
        set(gcf,'outerposition',[0 ss(4)/2 ss(3)/2 ss(4)/2])
        subplot(1,2,1); hold on;
        abs_df=abs(df);
        abs_df_smooth=abs(df_smooth);
        h1=plotplus(qnorm_aa/(4*pi), abs_df, adata0,'x');
        if sum(index_rejection)>0
            plotplus(qnorm_aa(index_rejection)/(4*pi),abs_df(index_rejection),adata1,'ok');
        end
        h2=plot(qnorm_aa(ind_sort)/(4*pi),abs_df_smooth(ind_sort),'-m');
        h3=plot(qnorm_aa(ind_sort)/(4*pi),df_mw(ind_sort),'-g');
        
        box on
        ylabel('abs(Deconvolution factor)','fontsize',14)
        xlabel('sin(\theta)/\lambda (Å^{-1})','fontsize',14)
        legend([h1 h2 h3], '\it{df}\rm{(\bf{H})}', 'smooth \it{df}', 'MW \it{df}','location', 'northwest')
        xlim([xmin,xmax])
        subplot(1,2,2); hold on
        phase_df=atan2(imag(df),real(df));
        phase_df_smooth=atan2(imag(df_smooth),real(df_smooth));
        h1=plotplus(qnorm_aa/(4*pi), phase_df, adata0,'x');
        if sum(index_rejection)>0
            plotplus(qnorm_aa(index_rejection)/(4*pi),phase_df(index_rejection),adata1,'ok');
        end
        h2=plot(qnorm_aa(ind_sort)/(4*pi),phase_df_smooth(ind_sort),'-m');
        ylabel('phase(Deconvolution factor)','fontsize',14)
        xlabel('sin(\theta)/\lambda (Å^{-1})','fontsize',14)
        xlim([xmin,xmax])
    end
end

if input_structure.verbose>0
    
    if input_structure.verbose >2
        f2=figure;
        max_deviation=max(abs(F_NDD_deviation));
        min_deviation=min(abs(F_NDD_deviation));
        if is_centrosymmetric
            hold on
            plotplus(qnorm_aa/(4*pi),abs(F_NDD_deviation),adata0,'x');
            if sum(index_rejection)>0
                plotplus(qnorm_aa(index_rejection)/(4*pi),F_NDD_deviation(index_rejection),adata1,'ok');
            end
            box on
            ylabel('F_{obs}^{NDD}/F_{calc}^{NDD}','fontsize',14)
            xlabel('sin(\theta)/\lambda (Å^{-1})','fontsize',14)
            set(gcf,'outerposition',[800    50   576   459])
            set(gcf,'outerposition',[ss(3)/2 ss(4)/2 ss(3)/2 ss(4)/2])
            
            line([min_stol max_stol],    repmat(input_structure.outlier_limit,1,2),'linestyle','--','color','m')
            line([min_stol max_stol], 1./repmat(input_structure.outlier_limit,1,2),'linestyle','--','color','m')
            text(min_stol,   input_structure.outlier_limit-0.2,     'Outlier limit')
            text(min_stol,1./input_structure.outlier_limit+0.2, '1/(Outlier limit)')
            legend('data','Outlier limit')
            xlim([xmin,xmax])
            if min_deviation < 1./input_structure.outlier_limit; ymin=min_deviation;     else ymin=0; end
            if max_deviation > input_structure.outlier_limit; ymax=max_deviation*1.05;     else ymax=input_structure.outlier_limit*1.05; end
            
            ylim([ymin, ymax])
            
        else
            
            subplot(1,2,1)
            
            hold on
            plotplus(qnorm_aa/(4*pi),abs(F_NDD_deviation),adata0,'x');
            if sum(index_rejection)>0
                plotplus(qnorm_aa(index_rejection)/(4*pi),abs(F_NDD_deviation(index_rejection)),adata1,'ok');
            end
            
            box on
            ylabel('abs(F_{obs}^{NDD}/F_{calc}^{NDD})','fontsize',14)
            xlabel('sin(\theta)/\lambda (Å^{-1})','fontsize',14)
            set(gcf,'outerposition',[ss(3) ss(4)/2 ss(3)/2 ss(4)/2])
            
            line([min_stol max_stol],    repmat(input_structure.outlier_limit,1,2),'linestyle','--','color','m')
            line([min_stol max_stol], 1./repmat(input_structure.outlier_limit,1,2),'linestyle','--','color','m')
            text(min_stol,   input_structure.outlier_limit-0.2,     'Outlier limit')
            text(min_stol,1./input_structure.outlier_limit+0.2, '1/(Outlier limit)')
            legend('data','Outlier limit')
            xlim([xmin,xmax])
            if min_deviation < 1./input_structure.outlier_limit; ymin=min_deviation;     else ymin=0; end
            if max_deviation > input_structure.outlier_limit; ymax=max_deviation*1.05;     else ymax=input_structure.outlier_limit*1.05; end
            ymin
            ymax
            ylim([ymin, ymax])
            subplot(1,2,2)
            hold on
            phase_deviation=atan2(imag(F_NDD_deviation),real(F_NDD_deviation));
            plotplus(qnorm_aa/(4*pi),phase_deviation,adata0,'x');
            if sum(index_rejection)>0
                plotplus(qnorm_aa(index_rejection)/(4*pi),phase_deviation(index_rejection),adata1,'ok');
            end
            box on
            ylabel('phase(F_{obs}^{NDD}/F_{calc}^{NDD})','fontsize',14)
            xlabel('sin(\theta)/\lambda (Å^{-1})','fontsize',14)
            set(gcf,'outerposition',[0 ss(4)/2 ss(3)/2 ss(4)/2])
            xlim([xmin,xmax])
        end
        
    end
    
end

sigma_obs_NDD=sigma_obs_NDD.*empiric_weight_correction;
A_obs_NDD=real(F_obs_NDD); B_obs_NDD=imag(F_obs_NDD);


%Avoid problems with infinite sigmas
ind_inf=isinf(sigma_obs_NDD);
sigma_obs_NDD(ind_inf)=max(sigma_obs_NDD(ind_inf~=1));
ind_nan=isnan(sigma_obs_NDD);
sigma_obs_NDD(ind_nan)=max(sigma_obs_NDD(ind_inf~=1));


if input_structure.verbose > 2
    
    
    adata0.hkl=hkl_baymem;
    adata0.sinlambda=qnorm_aa/(4*pi);
    
    f3=figure; hold on
    plotplus(qnorm_aa/(4*pi),A_calc_EDD,adata0,'x')
    plotplus(qnorm_aa/(4*pi),A_calc_NDD,adata0,'xr')
    plotplus(qnorm_aa/(4*pi),A_obs_EDD,adata0,'ob')
    plotplus(qnorm_aa/(4*pi),A_obs_NDD,adata0,'or')
    
    if sum(index_rejection(:))~=0
        adata1.hkl=hkl_baymem(index_rejection,:);
        adata1.sinlambda=qnorm_aa(index_rejection,:)/(4*pi);
        plotplus(qnorm_aa(index_rejection)/(4*pi),A_calc_EDD(index_rejection),adata1,'ok')
        plotplus(qnorm_aa(index_rejection)/(4*pi),A_calc_NDD(index_rejection),adata1,'ok')
        plotplus(qnorm_aa(index_rejection)/(4*pi),A_obs_EDD(index_rejection),adata1,'ok')
        plotplus(qnorm_aa(index_rejection)/(4*pi),A_obs_NDD(index_rejection),adata1,'ok')
    end
    
    output=[hkl_model qnorm_aa A_calc_EDD A_calc_NDD A_obs_EDD A_obs_NDD];
    fid=fopen('output.dat','w');
    fprintf(fid,'%4s%4s%4s%16s%16s%16s%16s%16s\n','h','k','l','qnorm_aa', 'A_calc_EDD', 'A_calc_NDD', 'A_obs_EDD', 'A_obs_NDD');
    fprintf(fid,'%4i%4i%4i%16.8e%16.8e%16.8e%16.8e%16.8e\n',output');
    fclose(fid);
    legend('EDD - calc','eNDD - calc','EDD - obs','eNDD - obs','location','northeastoutside')
    xlabel('sin(\theta)/\lambda (Å^{-1})','fontsize',14)
    ylabel('Structure factor','fontsize',14)
    xlim([xmin,xmax])
      set(gcf,'outerposition',[ss(3)/2 0 ss(3)/2 ss(4)/2])
    
end





if Uiso~=0
    disp('Deconvolving the the thermal motion')
    deconvl_type='_tm_deconvl';
    Umat=eye(3)*Uiso;
    deconvl_fac=DebyeWaller(M,Umat,hkl_baymem')';
    A_obs_NDD=A_obs_NDD./deconvl_fac;
    B_obs_NDD=B_obs_NDD./deconvl_fac;
    sigma_obs_NDD=sigma_obs_NDD./deconvl_fac;
    
end
number_of_electrons=A_obs_EDD(1);

if input_structure.verbose>0
    disp(['Assumes that the number of electrons is: ' num2str(number_of_electrons)])
end

%remove rejected reflections from lists
F_obs_NDD(index_rejection)=[];
F_calc_NDD(index_rejection)=[];
F_obs_EDD(index_rejection)=[];
F_calc_EDD(index_rejection)=[];
hkl_model(index_rejection,:)=[];

A_calc_EDD(index_rejection)=[];
A_calc_NDD(index_rejection)=[];
A_obs_EDD(index_rejection)=[];
A_obs_NDD(index_rejection)=[];
B_calc_EDD(index_rejection)=[];
B_calc_NDD(index_rejection)=[];
B_obs_EDD(index_rejection)=[];
B_obs_NDD(index_rejection)=[];
sigma_obs_NDD(index_rejection)=[];
hkl_baymem(index_rejection,:)=[];

% whos inputfile outputfile hkl_baymem A_obs_NDD B_obs_NDD sigma_obs_NDD
%Write to BayMEM file
writeFobs_BayMEM([inputfile],input_structure.outputfile,hkl_baymem,A_obs_NDD,B_obs_NDD,sigma_obs_NDD)
writeElectrons_BayMEM(input_structure.outputfile,input_structure.outputfile,number_of_electrons)
if input_structure.verbose>0
    disp(['Results written to: ' input_structure.outputfile])
    disp('')
    disp('--------------------------------------------------------')
    disp(sprintf('%.0f reflections were included',length(hkl_baymem)))
end

end

function does_function_exist(function_name)
if exist(function_name,'file')~=2
    error(['The function ' function_name ' was not found'])
end
end

function index_sigma_reject=reject_based_on_sigma(F_obs,sigma_obs,input_structure)
%3sigma criteria
sigma_Int=2*sigma_obs.*abs(F_obs);
index_sigma_reject=abs(F_obs).^2<=input_structure.sigma_reject*sigma_Int;
end

function print_sigma_rejected_reflections(hkl,F,sigma)
n_sigma=length(F);
if n_sigma>0
    fprintf('Found %i reflections below 3 sigma:\n',n_sigma)
    for i=1:n_sigma
        fprintf('hkl: %i %i %i  F_obs: %.3f , sigma: %.3f\n',hkl(i,:),F(i),sigma(i))
    end
end
end

function [F_obs_NDD,empiric_weight]=outliers_compared_to_model(F_obs_NDD,F_calc_NDD,input_structure)
F_NDD_deviation=F_obs_NDD./F_calc_NDD;
ind=F_NDD_deviation<1;
F_NDD_deviation(ind)=F_NDD_deviation(ind).^-1;
ind_outliers=F_NDD_deviation>input_structure.outlier_limit | F_NDD_deviation<1./input_structure.outlier_limit;

empiric_weight=F_NDD_deviation.^(input_structure.weight_exponent);

disp('--------------------------------------------------------')

switch input_structure.outlier_rule
    case 'none'
        fprintf('The outliers are kept. This may cause the MEM calculation to become ill behaved\n')
    case 'calc'
        F_obs_NDD(ind_outliers)=F_calc_NDD(ind_outliers);
        
    case 'reject'
        F_obs_NDD(ind_outliers)=NaN;
end
end


function print_outlier_reflections(hkl, deviation,input_structure)

n_outliers=length(deviation);
if n_outliers>0
    if input_structure.verbose > 1;
        fprintf('Found %i outliers listed below:\n',n_outliers);
    end
    
    if input_structure.verbose > 1
        for i=1:n_outliers
            disp(sprintf('hkl: %i %i %i Deviation factor: %f',hkl(i,:),deviation(i)))
        end
    end
    
    switch input_structure.outlier_rule
        case 'none'
            disp(sprintf('\nThe outliers are kept. This may cause the MEM calculation to become ill behaved'))
        case 'calc'
            disp(sprintf('\nF_calc_NDD are used instead of F_obs_NDD'))
        case 'reject'
            disp(sprintf('\nThese reflections are not used'))
        otherwise
            error('Unknown outlier rule given. Valid options are: ''none'', ''calc'', ''reject''')
    end
end
end

function [F_calc, A, B,ff,Z, hkl_model]=read_partial_F(filebase,element,M)
for k=1:length(element)
    atom_type=regexprep(element{k},'(\<[a-z])','${upper($1)}'); %capitalize first letter
    [a,  b]=get_f_coefficients(atom_type);
    %load elementwise structure factors
    filename=[filebase '_' element{k} '.m80'];
    try
        dat=importdata(filename);
    catch
        error(['Could not find file: ' filename ])
    end
    
    hkl_model=dat(:,1:3);
    
    F_calc(:,k)=dat(:,7);
    A(:,k)=dat(:,8);
    B(:,k)=dat(:,9);
    
    qnorm_aa(:,k)=sqrt(sum((M*hkl_model').^2,1))';
    ff(:,k)=sum(repmat(a,size(qnorm_aa(:,k),1),1).*exp(-(qnorm_aa(:,k)./(4*pi)).^2*b),2);
    Z(:,k)=sum(a,2);
    
end

stoichiometry=F_calc(1,:)./Z;
formula='';
for k=1:length(element)
    formula=[formula element{k} sprintf('%.2f',stoichiometry(k))];
end
disp(['The partial structure factors corresponds to the formula: ' formula])
end