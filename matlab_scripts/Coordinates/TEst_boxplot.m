clear
close all
clc

figure,

% colors = [1 0 0; 1 0 0; 0 0 1; 0 0.5 0; 0 0.5 0; 0 0.5 0];
colors = viridis(6);

x = boxplot(rand(100,6), 'Colors',colors);
% findall is used to find all the graphics objects with tag "box", i.e. the box plot
% hLegend = legend(flip(findall(gca,'Tag','Box')), {'Group A','Group B','Group C'})
all_items = findall(gca,'Tag','Box');
hLegend = legend([all_items(6);all_items(4);all_items(2)], {'Group A','Group B','Group C'})
