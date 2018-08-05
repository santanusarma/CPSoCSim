function view_thermal_matrices(A_HS, B_HS)

%Display thermal matrices
%Load thethermal matrices
figure(1)
set (1, 'color', [1 1 1])
subplot(2,2,1)
spy(A_HS)
subplot(2,2,2)
%spy(A_HZ)
heatmap(A_HS, [],[],'%0.2f', 'Colormap', 'money', 'Colorbar', true)
subplot(2,2,3)
spy(B_HS)
subplot(2,2,4)
heatmap(B_HS, [],[],'%0.2f', 'Colormap', 'money', 'Colorbar', true)


figure(2)
heatmap(A_HS, [],[],'%0.2f', 'Colormap', 'money', 'Colorbar', true)