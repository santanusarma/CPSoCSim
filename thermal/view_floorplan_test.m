% Floor Plan Viewing  Test Cases (Pass)
% function view_floorplan(w, h, x, y, hFigure, label, z)

%==========================================================================
%Version    : 1.0
%Date       : July 2012
%Author     : Santanu Sarma, University of California, Irvine
%Address    : 3069, Embedded System Lab, Bren Hall, UCI, CA 92697
%Email      : santanus@uci.edu
%Web        : https://students.ics.uci.edu/~santanus/
%==========================================================================

 
 clear all; clc
 
 
 % Load the EV6 fllor plan 
[unit_name_1, width_1, height_1, x_left_1, y_bottom_1] = ...
    textread('ev6_1x1.flp', '%s %f %f %f %f');
 
% Inputs:
w=width_1
h=height_1
x=x_left_1
y=y_bottom_1
hFigure=1
label=unit_name_1
z= randi(256, length(width_1),1 )

% Floor plan with out lables and colors (Pass)
view_floorplan(w, h, x, y, hFigure)


% Floor plan with lables  (Pass)
view_floorplan(w, h, x, y, 2, label)


% Floor plan withlables and colors (Pass)
view_floorplan(w, h, x, y, 3, label,z)


%% Floor plan with lables and gray colors for multi cores (Pass)
[label_mc, w_mc, h_mc, x_mc, y_mc] = ...
    textread('mc_ev6_2x2.flp', '%s %f %f %f %f');
z_mc= randi(256, length(w_mc),1 );
view_floorplan(w_mc, h_mc, x_mc, y_mc, 4, label_mc,z_mc)



%floor plan without lables for multicores (Pass)
nBlocks=length(label_mc)
for i=1:nBlocks
label_mc{i}='';
end
view_floorplan(w_mc, h_mc, x_mc, y_mc, 5, label_mc,z_mc)



%% Floor plan with lables and gray colors for 16 cores multi cores (Pass)
cd models/4x4/
[label_mc, w_mc, h_mc, x_mc, y_mc] = ...
    textread('ev6_4x4.flp', '%s %f %f %f %f');
z_mc= randi(256, length(w_mc),1 );
for i=1:length(label_mc)
label_mc{i}='';
end
view_floorplan(w_mc, h_mc, x_mc, y_mc, 5, label_mc,z_mc)


%floor plan with colored floor plan (RGB) % TBD



%==========================================================================
clear all; close all;
% Hotspot Floor Plan of EV6 with 30 blocks
%               w   h   x   y
FloorPlanEv6=[0.0049	0.0062	0	0.0098
0.016	0.0098	0	0
0.0049	0.0062	0.0111	0.0098
0.0031	0.0026	0.0049	0.0098
0.0031	0.0026	0.008	0.0098
0.001033	0.0007	0.0049	0.0124
0.001033	0.0007	0.005933	0.0124
0.001033	0.0007	0.006967	0.0124
0.001033	0.0007	0.008	0.0124
0.001033	0.0007	0.009033	0.0124
0.001033	0.0007	0.010067	0.0124
0.0011	0.0009	0.0049	0.0131
0.0011	0.0009	0.006	0.0131
0.00055	0.00038	0.0049	0.014
0.00055	0.00038	0.00545	0.014
0.00055	0.00038	0.006	0.014
0.00055	0.00038	0.00655	0.014
0.0011	0.00095	0.0049	0.01438
0.0011	0.00095	0.006	0.01438
0.0011	0.00067	0.0049	0.01533
0.0011	0.00067	0.006	0.01533
0.0009	0.00135	0.0071	0.01465
0.0013	0.00135	0.008	0.01465
0.0009	0.00067	0.0093	0.01533
0.0009	0.00067	0.0102	0.01533
0.0018	0.00223	0.0093	0.0131
0.0009	0.00155	0.0071	0.0131
0.0013	0.00095	0.008	0.0137
0.00065	0.0006	0.00865	0.0131
0.00065	0.0006	0.008	0.0131]

w=FloorPlanEv6(:,1);
h=FloorPlanEv6(:,2);
x=FloorPlanEv6(:,3);
y=FloorPlanEv6(:,4);

%FP=[w h x y]

%% Calculate total width and height
x_left = min(x);
x_right = max(x + w);
W = x_right - x_left;
y_bottom = min(y);
y_top = max(y + h);
H = y_top - y_bottom;

%Temeperature Steady State  form HotSpot
TempTrace=[60.09	60.08	60.09	61.72	63.26	63.2	63.22	63.16 ...
    60.31	60.27	60.28	61.15	61.14	61.49	61.51	61.51 ...	
    61.45	61.07	61.08	60.11	60.13	61.42	60.5	67.16 ...
    67.25	63.23	60.2	64.7	61.05	60.97]

%       Die     TIM     HSep    HS      
TempSS=[51.2	51.14	50.95	50.82
50.71	50.66	50.49	50.41
51.53	51.47	51.26	51.07
57.94	56.93	53.58	52.16
62.64	60.81	54.69	52.5
60.62	58.95	53.41	51.92
61.05	59.35	53.69	52.09
60.87	59.26	53.89	52.29
55.65	55.24	53.88	52.41
55.53	55.16	53.94	52.43
55.21	54.85	53.63	52.23
56.12	55.38	52.91	51.8
56.4	55.67	53.23	52.02
55.96	55.17	52.55	51.64
56.17	55.36	52.68	51.71
56.39	55.58	52.9	51.86
56.3	55.54	53	51.93
54.82	54.22	52.25	51.52
55.24	54.62	52.59	51.73
52.32	52.17	51.69	51.28
52.86	52.66	52.02	51.46
56.14	55.36	52.76	51.78
55.36	54.85	53.16	51.91
70.58	66.76	54.05	51.69
71.14	67.2	54.05	51.6
62.3	60.46	54.31	52.17
54.71	54.39	53.31	52.14
64.35	62.01	54.2	52.24
57.53	56.71	53.98	52.33
56.92	56.21	53.82	52.31]



%Make an Image from the block information i.e 1-D immage coordinates to 
%2-D image coordinates

%Step 0: Normalize the image floor plan; find the min dim of hight and width 
%Step 2: make sure that the size of the image is multiple of each block
%Step 3: make equivalent rectangle image for each block with the given
%temperature of the block
%Step 4: convert each block color as per the temperature
%Step 5: 

min_dim=min(min(FloorPlanEv6(:,1:2)))

FloorPlanEv6Norm=FloorPlanEv6./min_dim

%Step 1: define the pixle size and dimension of the image
Wp=1024
Hp=1024
N=Wp*Hp  %No of elements in the thermal map 
%Wp, Hp are no of pixles that has to map to the no 

%precision=min_dim/10;
precision=max(W,H)/max(Wp,Hp)
i=1
NoIndx_x=round(w(i)/precision)
NoIndx_y=round(h(i)/precision)

FloorPlanEv6Indx=ceil(FloorPlanEv6./precision)
wx=FloorPlanEv6Indx(:,1);
hx=FloorPlanEv6Indx(:,2);
xx=FloorPlanEv6Indx(:,3);
yx=FloorPlanEv6Indx(:,4);

% NoIndx_x=round(W/precision)
% NoIndx_y=round(H/precision)

I=zeros(Wp,Hp);
for bn=1:30  % for each block
    
    for i=1:wx(bn)+1  % overlap one pixel 
        for j=1: hx(bn)+1
            
        Indx_x= (xx(bn)+ i);
        Indx_y= (yx(bn)+ j);
        if (Indx_x<=Wp)&(Indx_y<=Hp)
        I(Indx_x,Indx_y)=TempSS(bn,1);
        end
        end
    end
end


I1=imrotate(I,90); %rotate the image  by 90 degree
size(I1)
[B lims map] = real2rgb(I1, 'jet1');
hIm = imshow(B);
set(gcf, 'Colormap', map);
set(gca, 'CLim', lims);
set(hIm, 'CDataMapping', 'scaled');
colorbar;

%--------------------------------------------------------------------------
% Check thermal map with spectral coeficient truncation 

Tss=TempSS(:,1) %Steady state temperature 
NoOfLoc=length(Tss)
Tss_dctcoef=dct(Tss)
NoOfCoef=length(Tss_dctcoef)
D=dctmtx(NoOfLoc)
Dinv=inv(D);
Tss_dctcoef1=D*Tss

Tss_hat=idct(Tss_dctcoef)
Tss_hat1=inv(D)*Tss_dctcoef1 % Tss_hat=Dinv*Tss_dctcoef1
Tss_Err=Tss_hat-Tss_hat1
Tss_CoefErr=Tss_dctcoef-Tss_dctcoef1


[Tss_dct1a, I1]=sort(Tss_dctcoef1,'descend')

[Tss_dct1b, I2]=sort(abs(Tss_dctcoef1),'descend')
norm1a=sum(abs(Tss_dctcoef1)) % energy ??
norm1b=norm(Tss_dctcoef1,1)
energy=sum(Tss_dctcoef1.^2)
n=20
%SensLoc=randi(30,5,1);
SensLoc=I2(1:n);
for n=1:NoOfCoef; %no of coefficient
CoefLoc=I2(1:n)
%SensLoc=CoefLoc;

%SensLoc=randi(30,n,1);
Tss_dctcoef_truc=Tss_dctcoef1(CoefLoc)
Phi=Dinv(SensLoc,CoefLoc)
cond_no(n)=cond(Phi)
CLSE=inv(Phi'*Phi)*Phi'*Tss(SensLoc)
Tss_hat=Dinv(:,CoefLoc)*Tss_dctcoef_truc
Tss_hat1=Dinv(:,CoefLoc)*CLSE

Tss_err=Tss-Tss_hat;
Tss_err1=Tss-Tss_hat1;
MSE(n)=sum(Tss_err.^2)/(length(Tss_err)-1)
MSE1(n)=sum(Tss_err1.^2)/(length(Tss_err1)-1)

view_thermalmap('ev6_hotspot.flp',Tss_hat, 1024, 1024,n)
pause
end


figure(1)
plot([1:NoOfCoef],MSE), xlabel('No of Sensors')
ylabel('MSE'), axis tight, grid 
figure(2)
plot([1:NoOfCoef], MSE1)
figure(3)
plot([1:NoOfCoef], MSE,[1:NoOfCoef], MSE1 )

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
%STOMP
N=length(Tss)
C=zeros(N,1)
s=[] % coefficient location 
Threshold=50
n=4; % no of sensor
y=1:n;  % sensor location 
t=Tss(y)
r=t
Phi=inv(D);
iter=1
MaxIter=10000;
MSEtol=1e-5;
MSEpre=inf;
while iter<MaxIter
    rand_loc=randperm(N);
    y=rand_loc(1:n) 
   
    
    
    Cr=Phi(y,:)'*r
    [Cs, Ps]=sort(Cr,'descend')
    
    p=find(abs(Cr)>(Threshold))
    s1=s
    s=[s; p]
    
    norm_r=norm(r)
    rank_phi=rank(Phi(y,s)'*Phi(y,s))
    cond_no=cond(Phi(y,s))

    if (rank_phi==length(s))
    C(s)=inv(Phi(y,s)'*Phi(y,s))*Phi(y,s)'*t
    Tss_hat=Phi(y,:)*C;
    r=t-Tss_hat;
    
    else % reject s
        %s=s(1:end-length(p))
        s=s1;
        
    end
    
    NoOfCoef=length(s)
    s
    y
    MSE=sum(r.^2)/length(r)
%     if (MSE<MSEpre)
%         MSEpre=MSE
%         ypre=y;
%     else
%         MSE=MSEpre
%         y=ypre
%         %s=s1
%     end
    iter=iter+1
    if MSE<MSEtol
        break
    end
    %pause
end

%--------------------------------------------------------------------------

SenLoc=I2(1:5)
Tss_m=Tss(SenseLoc)
Phi=D(SenLoc,:)
Phi_p=inv(Phi'*Phi)*Phi'
Phi_pinv=pinv(Phi)

Tss_coeff_lse=Phi_p*Tss_m

Tss_hat=idct(Tss_coeff_lse)

Tss_hat=pinv(Phi)*Tss_dct1b



%--------------------------------------------------------------------------
%Generate analytical thermal maps
clear all; close all
Wp=256
mag=0*sin(2*pi*0.01*[1:Wp]).*cos(2*pi*0.03*[1:Wp]);
Mag=repmat(mag, Wp,1);
size(Mag)
TM=2*peaks(Wp)+50;
TM1=TM+Mag;
[r, c]=size(TM1)
[B lims map] = real2rgb(TM1, 'jet2');
hIm = imshow(B);
set(gcf, 'Colormap', map);
set(gca, 'CLim', lims);
set(hIm, 'CDataMapping', 'scaled');
colorbar;

%Step 1: convert the image to 1-d vector 
[r c]=size(TM1)
for i=1:r
    for j=1:c
        indx=sub2ind([r c], i,j);
        TM_1d(indx)=TM1(i,j);

    end
end

%--------------------------------------------------------------------------
% k-LSE /dCT based - Sharif Reda
% http://dl.acm.org/citation.cfm?id=1837291



TM_1d=randi(100,16,1)
%TM_1d=magic(4)
n=length(TM_1d)
D=dctmtx(n)
%compute the dct 
TM_dct1=D*TM_1d;

TM_dct2=D*TM_1d*D'; %Some problem with dimension 

dct_err=TM_dct1-TM_dct2;


%--------------------------------------------------------------------------
% Eigen Maps Method of Image Reconstrcutions 


figure(2), plot(TM_1d)
%Step 2: compute the covarience of the TM_1d vector 
C1=cov(TM_1d)
C=cov(TM1);
size(C)

%Step 3: find the eigen values and eigen vectors
[V, D]=eig(C)

%Step 4: 

