% An example of how to use MagInset.m

clear all;

% **************************************
% * Example 1: Single axes on a figure *
% **************************************

% Create a figure with a sine wave with random noise on it:
h1 = figure;
title('Single Axes MagInset Example');
hold on;
f0 = 10000;
t = 0:(1/f0)/1000:((1/f0)*2);
y = 5*sin(2*pi*f0*t) + rand(size(t));
plot(t*1000, y);
grid on;
xlabel('Time (ms)');
ylabel('Amplitude (V)');
% It is important to set the figure size and limits BEFORE running
% MagInset!
xlim([t(1) t(end)]*1000);
ylim([-6 10]);
set(h1, 'Position',[105   647   560   420]);


% Once happy with your figure, add an inset:
MagInset(h1, -1, [0.04 0.045 1 4], [0.05 0.09 5.5 9], {'NW','NW';'SE','SE'});

% ****************************************
% * Example 2: Multiple axes in a figure *
% ****************************************

% Create a figure with a sine and cosine wave with random noise on it:
h2 = figure;

% Sine Wave:
ax(1) = subplot(2,1,1);
title(sprintf('MagInset Example with Subplots\n'));
hold on;
f0 = 10000;
t = 0:(1/f0)/1000:((1/f0)*2);
y = 5*sin(2*pi*f0*t) + rand(size(t));
plot(t*1000, y);
grid on;
xlabel('Time (ms)');
ylabel('Amplitude (V)');
xlim([t(1) t(end)]*1000);
ylim([-6 10]);

% Cosine Wave:
ax(2) = subplot(2,1,2);
hold on;
y = 5*cos(2*pi*f0*t) + rand(size(t));
plot(t*1000, y,'r');
grid on;
xlabel('Time (ms)');
ylabel('Amplitude (V)');
xlim([t(1) t(end)]*1000);
ylim([-6 10]);

% Link the x axes:
linkaxes(ax,'x');

% Set the figure size:
set(h2, 'Position',[682   381   560   685]);

% Once happy with your figure, add an inset:
% Cosine inset:
MagInset(h2, ax(1), [0.04 0.045 1 4], [0.05 0.09 5.5 9], {'NW','NW';'SE','SE'});
% Sine insets:
MagInset(h2, ax(2), [0.015 0.02 1.5 4], [0.038 0.078 5.5 9], {'NW','NW';'SE','SE'});
MagInset(h2, ax(2), [0.18 0.185 1.5 4], [0.125 0.165 5.5 9], {'NE','NE';'SW','SW'});