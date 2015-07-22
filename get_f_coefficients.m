function [a,b] = get_f_coefficients(atom_type,database)
% atom_type=regexprep(atom_type,'(\<[a-z])','${upper($1)}'); %capitalize first letter
if nargin < 2
    database='jana';
end

switch database
    case 'prior'
        try
            string= fileread('factors_atomsions.fit');
        catch err
            error('The file ''factors_atomsions.fit'' must be in Matlab''s search path, for get_f_coefficients to work')
        end
        temp=regexp(string,['[ \t]+' atom_type '[ \n\r]+\s+' repmat('[ \t]*([0-9\.]+)',1,6) '[ \n\r]+' repmat('[ \t]*([0-9\.]+)',1,6)],'tokens');
        ffc=str2double(temp{1});
        
        a=ffc(1:6); b=ffc(7:12);
    case 'jana'
        try
            load('jana_atomic_database.mat','database')
        catch err
            error('The file ''jana_atomic_database.mat'' must be in Matlab''s search path, for get_f_coefficients to work')
        end
        ffc=database.(atom_type).xray_form_factor_analytical;
        a=[ffc(1:2:9)];
        b=[ffc(2:2:8) 0];
end
end
