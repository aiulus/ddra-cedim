function s = num2str_cell(x)
    if ischar(x) || isstring(x), s = char(x);
    elseif islogical(x), s = char(string(x));
    else, s = num2str(x, '%.10g'); end
end