clear all
close all

%%

Ncom = 100;
Nrad = 100;
t = 0:(Nrad-1);
xrad = exp(1j*pi*t.^2/Nrad);

Fr=gallery('circul',[xrad.'; zeros(Ncom-Nrad,1)]);

imax = 15;
jmax = 20;
iter = 500;

tic

imp = zeros(jmax,imax+1);
avgpar = zeros(jmax,imax+1);
avgevm = zeros(jmax,imax+1);

for ii = 0:imax
   
    if ii == 0
        S = [];
    else
        S = 1:round(Ncom/ii):Ncom;
    end
    if length(S) ~= ii
        S(end) = [];
    end
    St = 1:Ncom;
    St(S) = [];
    
    for jj = 1:jmax

        par = zeros(iter,1);
        newpar = zeros(iter,1);
        evm = zeros(iter,1);
        for k = 1:iter
            xcom = randsrc(Ncom,1,[1 -1 1j -1j]);
            one = xcom==1;one(S)=0;
            negone = xcom==-1;negone(S)=0;
            onej = xcom==1j;onej(S)=0;
            negj = xcom==-1j;negj(S)=0;
            
            cvx_begin quiet
                cvx_precision(1e-1)
                variable x(Ncom) complex
                minimize( norm(Fr*(xcom+x),Inf) )
                subject to
%                     norm(x(S))/norm(xcom(S)) <= 1000;
%                     norm(x(St))/norm(xcom(St)) <= sst+ssp*jj;
                    abs(sum(x(one))) <= 0.01;
                    abs(sum(x(negone))) <= 0.01;
                    abs(sum(x(onej))) <= 0.01;
                    abs(sum(x(negj))) <= 0.01;
                    
                    %Circular 
                    abs(x(St)) <= .0175*jj;

%                     %square                
%                     abs(real(x)) <= .05+.01*jj;
%                     abs(imag(x)) <= .05+.01*jj;

            cvx_end

            x2 = (xcom+x);

            zcom = Fr*xcom;
            zcom2 = Fr*x2;

            par(k) = 10*log10(norm(zcom,Inf)^2/(norm(zcom)^2/Ncom));
            newpar(k) = 10*log10(norm(zcom2,Inf)^2/(norm(zcom2)^2/Ncom));
            evm(k) = norm(x(St))/norm(xcom(St));
        end
        avgevm(jj,ii+1) = mean(evm);
        imp(jj,ii+1) = mean(par) - mean(newpar);
        avgpar(jj,ii+1) = mean(newpar);
        stdpar(jj,ii+1) = std(newpar);
        fprintf('time: %g iter: %g\n',toc,(jmax*ii+jj)/((jmax)*(imax+1))*100);
    end
    
end
%%
xind = 0:imax; yind =mean(avgevm,2)';
save(strcat('C:\MatlabData\Data',datestr(now,30)),'imp','avgpar','*ind','stdpar');

%%
xind_i = interp1(xind,1:.01:length(xind),'linear');
yind_i = interp1(yind,1:.01:length(yind),'linear');

[gridXi, gridYi] = meshgrid(xind_i,yind_i);
[gridX, gridY] = meshgrid(xind,yind);

mkedB_i = interp2(gridX,gridY,avgpar,gridXi,gridYi,'spline');

figure (350); set(gcf,'units','inches');
set(gca,'FontSize',20); set(findall(gcf,'type','text'),'FontSize',20);
set(gcf,'Position',[2.5833    3.9271    3.3750    4.0729]);
imagesc( xind_i, yind_i, mkedB_i); hold on;
caxis([1 7]);cc=colorbar;
set(gca,'FontSize',20)
title('Constraints for Radar Spreading PAR Reduction (Circular Boundary)')
ylabel(cc,'PAR (dB)')
% surf(gridXi,gridYi,RMSEdB_i); colorbar;
set(gca,'YDir','normal'); %flip the y axis, default for imagesc is 'reverse'
[Xg, Yg] = ndgrid(xind,yind);
plot(Xg,Yg,'kx'); %superimpose the uninterp'd data
xlabel('Number of Large Constraint Symbols'); ylabel('Max EVM of Small Constraint Symbols'); 
grid on
hold off
