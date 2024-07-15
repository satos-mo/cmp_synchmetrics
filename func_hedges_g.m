function [d,g] = func_hedges_g(x, y)
% [Cohen_d,Hedge_g] = hedges_g(sample1, sample2);

    mean_x = mean(x,'omitnan');
    mean_y = mean(y,'omitnan');
    
    var_x = var(x,'omitnan');
    var_y = var(y,'omitnan');
    
    n_x = sum(~isnan(x));
    n_y = sum(~isnan(y));
    
    pooled_std = sqrt(((n_x-1)*var_x+(n_y-1)*var_y)/(n_x+n_y-2));    
    d = (mean_x - mean_y) / pooled_std;    
    correction_factor = 1-(3/(4*(n_x+n_y)-9));
    g = d * correction_factor;
end


