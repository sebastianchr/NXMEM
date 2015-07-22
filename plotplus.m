% First plot the data
function hLine=plotplus(varargin)

x=varargin{1};
y=varargin{2};
adata=varargin{3};
if nargin>3
    hLine = plot(x, y,varargin{4:end});
else
    hLine = plot(x, y);
end

setappdata(hLine,'additional_data',adata)

dcm = datacursormode(gcf);
datacursormode on
set(dcm,'updatefcn',@myfunction)
%
% hTarget = handle(hLine);
% hDatatip = dcm.createDatatip(hTarget);

end

function output_txt = myfunction(obj,event_obj)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% add_data     Structure containing addtional information for DataTip
% output_txt   Data cursor text string (string or cell array of strings).
pos = get(event_obj,'Position');
ind = get(event_obj,'DataIndex');

temp = {['X: ',num2str(pos(1),4)],...
    ['Y: ',num2str(pos(2),4)]};
%loop through additional data


try
    p=get(event_obj,'Target');
    add_data=getappdata(p,'additional_data');
    if  isfield(add_data,'format');
        format_struct=add_data.format;
        add_data=rmfield(add_data,'format');
        isformat=1;
    else
        isformat=0;
    end
    fields = fieldnames(add_data);
    for i = 1:numel(fields)
        val=add_data.(fields{i})(ind,:);
        format='%.2f';
        if isformat
            if isfield(format_struct,fields{i})
                format=format_struct.(fields{i});
            end
        end
        formatstr=repmat([format ', '],1,numel(val));
        if iscell(val)
        str=[fields{i} ': ' sprintf(formatstr(:),val{:})];
        else
            str=[fields{i} ': ' sprintf(formatstr(:),val(:))];
        end
        temp{end+1}=str;
        
    end
    output_txt=temp;
catch ME
    disp(ME)
end
end

