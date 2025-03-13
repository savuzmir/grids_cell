function out = print_ttest(a, varargin)

if isempty(varargin)
    [~, b, ~, d] = ttest(a);
else
    if varargin{1} == 1
        [~,b,~,d] = ttest(a, varargin{2});
    else
        [~,b,~,d] = ttest2(a, varargin{2});
    end
end

out = [b, d.tstat, d.df];

end