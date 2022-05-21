close all
clc
clear
figure,
x = [1,2,3,4,5,1,2,3,4,6];
group = [1,1,2,2,2,3,3,3,4,4];
positions = [1 1.25 2 2.25];
boxplot(x,group, 'positions', positions);

set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) ])
set(gca,'xticklabel',{'Direct care','Housekeeping'})

color = ['c', 'y', 'c', 'y'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:2), 'Feature1', 'Feature2' );