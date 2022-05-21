function [] = save_PDF(fig_name,pdf_Name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    DPI = 600; % Dots per square inch.  Higher dpi will give higher resolution
    set(0, 'CurrentFigure', fig_name)
    %This will save a "cropped" pdf in the directory
    set(gcf,'Units','inches')
    h=get(gcf,'Position');
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0, 0 ,h(3), h(4)]);
    set(gcf, 'PaperSize', [h(3), h(4)])
    print('-dpdf',strcat('-r',num2str(DPI)),pdf_Name)
end

