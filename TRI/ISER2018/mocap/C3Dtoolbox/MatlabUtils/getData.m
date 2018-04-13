% ========================================================================
% getData(label, out)
% ========================================================================
% Gets the data belonging to a particular label from the extracted dataset

function data = getData(label, out)

l = length(out.labels);
for i = 1:l
    if strcmpi(out.labels{i}, label)
        data = out.data(:,i);
        return
    end
end
data = 0;