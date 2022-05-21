function [] = cdfplot2(error)
index_error = isnan(error);
error(index_error) = Inf;
cdfplot(error(:));
end