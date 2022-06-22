function h = cdfplot2(error)
index_error = isnan(error);
error(index_error) = Inf;
h = cdfplot(error(:));
end