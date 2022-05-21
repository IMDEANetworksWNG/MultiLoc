close all
x = [1:100];
y = [1:100];

std_noisy = 1.2;
x_noisy = x + randn(size(x))*std_noisy;
y_noisy = y + randn(size(y))*std_noisy;

figure
plot(x, y, 'gd');
hold on
plot(x_noisy, y_noisy)

error_noisy = sqrt((x_noisy-x).^2+(y_noisy-y).^2);

%figure
%cdfplot(error_noisy(:))

x_k = x_noisy(1);
y_k = y_noisy(1);
initialState = [x_k;0;y_k;0];
KF = trackingKF('MotionModel','2D Constant Velocity','State',initialState);

T  = 1.5;

pos = [x_noisy', y_noisy'];

figure
hold on
for k = 1:size(x_noisy,2)
    
    pstates(k,:) = predict(KF, T);
    cstates(k,:) = correct(KF, pos(k,:));
end

plot(pos(:,1),pos(:,2),'-k.', pstates(:,1),pstates(:,3),'-+', ...
    cstates(:,1),cstates(:,3),'-o')
xlabel('x [m]')
ylabel('y [m]')
grid
xt  = [x-2 pos(1,1)+0.1 pos(end,1)+0.1];
yt = [y pos(1,2) pos(end,2)];
hold on
plot(x, y, 'gd');

legend('Object position', 'Corrected position', 'GT')
%text(xt,yt,{'First measurement','First position','Last position'})

error_corr = sqrt((cstates(:,1)'-x).^2+(cstates(:,3)'-y).^2);
error_p = sqrt((pstates(:,1)'-x).^2+(pstates(:,3)'-y).^2);
figure
cdfplot(error_corr(:))
hold on
cdfplot(error_noisy(:))
cdfplot(error_p(:))
legend("Corr", "Noisy", "predicted")
