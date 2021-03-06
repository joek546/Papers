load ToneReservationVsEVM
% Figures for radar spreading PAR reduction paper
xind_i = interp1(xind,1:.01:length(xind),'linear');
yind_i = interp1(yind,1:.01:length(yind),'linear');

[gridXi, gridYi] = meshgrid(xind_i,yind_i);
[gridX, gridY] = meshgrid(xind,yind);

mkedB_i = interp2(gridX,gridY,avgpar,gridXi,gridYi,'spline');

Fig1 = figure (1); set(gcf,'units','inches');
set(gca,'FontSize',8); set(findall(gcf,'type','text'),'FontSize',8);
set(gcf,'Position',[22    5    3.375    3.375]);
Plot1 = imagesc( xind_i, yind_i, mkedB_i); hold on;
set(gca,'units','inches','Position',[0.5    0.4    2.4    2.8]);
% set(Plot1, 'Position', [0.5 0.38 2.5 3.0]);
caxis([1 7]);cc=colorbar;
set(cc,'units','inches', 'Position',[3, 0.4, 0.08, 2.8],'FontSize',8)
set(gca,'FontSize',8)
% title('Tone Reservation Vs. EVM Constraints (Cir. Boundary)')
ylabel(cc,'PAR (dB)')
% surf(gridXi,gridYi,RMSEdB_i); colorbar;
set(gca,'YDir','normal'); %flip the y axis, default for imagesc is 'reverse'
[Xg, Yg] = ndgrid(xind,yind);
plot(Xg,Yg,'kx','MarkerSize',3); %superimpose the uninterp'd data
xlabel('Number of Tone Reservation'); ylabel('Max EVM Constraint'); 
grid on;
saveas(Fig1,'TRvsEVM','fig');
matlabfrag('TRvsEVM')
hold off

% figure (2); set(gcf,'units','inches');
% set(gca,'FontSize',12); set(findall(gcf,'type','text'),'FontSize',12);
% set(gcf,'Position',[2.5833    3.9271    3.3750    4.0729]);
% imagesc( xind_i, yind_i, mkedB_i); hold on;
% caxis([1 7]);cc=colorbar;
% set(gca,'FontSize',12)
% title('Constraints for Radar Spreading PAR Reduction (Circular Boundary)')
% ylabel(cc,'PAR (dB)')
% % surf(gridXi,gridYi,RMSEdB_i); colorbar;
% set(gca,'YDir','normal'); %flip the y axis, default for imagesc is 'reverse'
% [Xg, Yg] = ndgrid(xind,yind);
% plot(Xg,Yg,'kx'); %superimpose the uninterp'd data
% xlabel('Number of Large Constraint Symbols'); ylabel('Max EVM of Small Constraint Symbols'); 
% grid on
% % saveas(gcf,'RMSE_vsN','fig');
% % matlabfrag('TRvsEVM')
% hold off