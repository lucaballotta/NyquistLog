function NyquistLog(sys,NgPlot,tol)

% NYQUISTLOG makes a polar plot of the open loop transfer function sys: the
% amplitude M=|sys| is converted to a logarithmic scale according to the
% following function:
%         _
%        |   M^(log10(2))       if M < 1
% L(M)=|
%        |_  2-M^(-log10(2))    if M > 1
%
% The optional boolean parameter NGPLOT specifies if also the diagram for
% negative frequencies should be plotted.
% The optional parameter TOL specifies the tolerance for quasi zero-pole
% cancelation by the MINREAL function.
% Default values are NGPLOT = TRUE and TOL = SQRT(EPS).

% May 2021
% Luca Ballotta
% Department of Information Engineering,
% University of Padova,
% Padova, Italy.
%
% Modified version of the function Closed_Logarithmic_Nyquist by
% Roberto Zanasi and Federica Grossi.
%
% May be distributed freely for non-commercial use,
% but please leave the above info unchanged, for
% credit and feedback purposes.

%%%%%%%%%
% Input check
% Checking system variable:
[hlp,den]=tfdata(sys,'v');
k=1;
while (hlp(k) == 0) k=k+1; end
num=(hlp(k:end));
hlp=size(size(den));
if (hlp(2)> 2)
    error('NyquistlLog:InputError',...
        'Only monovariable systems allowed.');
end
% Checking delay:
dly=get(sys,'ioDelay');
if (dly>0)
    error('NyquistlLog:InputError',...
        'Delay not allowed.');
end
% Checking continuous-time:
sorz=get(sys,'variable');
if (sorz ~= 's')
    error('NyquistlLog:InputError',...
        'Only continuous systems allowed.');
end
% Checking system order:
k=1;
while (den(k) == 0) k=k+1; end
den=(den(k:end));
sysdim=size(den); sysdim=sysdim(2)-1;
if (sysdim == 0)
    error('NyquistlLog:InputError',...
        'Denominator order of zero not allowed.');
end
numdim=size(num); numdim = numdim(2)-1;
if (numdim > sysdim)
    error('NyquistlLog:InputError',...
        'Only causal systems allowed.');
end
% Checking optional arguments
if nargin<2
    NgPlot = true;
elseif (NgPlot~=true && NgPlot~=false)
    error('NyquistlLog:InputError',...
        '''NgPlot'' must be boolean valued.');
end
if nargin<3
    tol = sqrt(eps);
elseif (not(isscalar(tol)) && tol<=0)
    error('NyquistlLog:InputError',...
        '''tol'' must be positive valued.');
end
% Checking open-loop stability
sys=minreal(sys,tol);
sys_poles = pole(sys);
if (any(real(sys_poles)>0))
    error('NyquistlLog:InputError',...
        'Only stable systems allowed.');
end

%%%%%%%%%
% Plot params
n=2;             % Base of the L function
Lw=1.2;          % Linewidth
Nd=4;            % Number of level lines
clear XPlot
XPlot.c='r'; XPlot.Lw=Lw; XPlot.NgPlot=NgPlot; XPlot.sys = sys;

%%%%%%%%%%%%%%%%%%%%%%%%
figure(); clf
plot(2*exp(1i*(0:0.01:1)*2*pi),'-')
hold on
grid on
plot(exp(1i*(0:0.01:1)*2*pi),'-')
plot(exp(1i*pi),'*r')
text(-1.05,-0.05,'-1')
xx=1*exp(1i*pi/4);
text(real(xx),imag(xx),[num2str(0) 'dB'])
for ii=1:Nd
    plot(1/(n^ii)*exp(1i*(0:0.01:1)*2*pi),':')
    plot((2-1/(n^ii))*exp(1i*(0:0.01:1)*2*pi),':')
    if ii<3
        xx=1/(n^ii)*exp(1i*pi/4);
        text(real(xx),imag(xx),[num2str(-ii*20) 'dB'])
        xx=(2-n^(-ii))*exp(1i*pi/4);
        text(real(xx),imag(xx),[num2str(ii*20) 'dB'])
    end
end
for ii=1:12 % plot of sectors of pi/6
    ps=2*exp(1i*ii*pi/6);
    plot([0 real(ps)],[0 imag(ps)],':')
end
plot([-1 1]*3,[-1 1]*0,':')
plot([-1 1]*0,[-1 1]*3,':')
axis equal
axis([-1 1 -1 1]*2.2)

[mult,freqs_inf] = groupcounts(imag(sys_poles(real(sys_poles)==0 & imag(sys_poles)>0)));
n_freqs_inf = length(freqs_inf);
if isempty(freqs_inf)
    [RE,IM] = nyquist(XPlot.sys,logspace(-Nd,Nd,1000));
else
    n_points = 1000/(n_freqs_inf+1);
    freq_inf_th = .01;
    while max(abs(freqresp(sys,freqs_inf-freq_inf_th))) <= 1
        freq_inf_th = freq_inf_th / 2;
    end
    [RE,IM] = nyquist(XPlot.sys,logspace(-Nd,log10(freqs_inf(1)-freq_inf_th),n_points));
    for ii = 1:n_freqs_inf-1
        [RE_i,IM_i] = nyquist(XPlot.sys,logspace(log10(freqs_inf(ii)+freq_inf_th),log10(freqs_inf(ii+1)-freq_inf_th),n_points));
        RE = cat(3,RE,RE_i);
        IM = cat(3,IM,IM_i);
    end
    [RE_f,IM_f] = nyquist(XPlot.sys,logspace(log10(freqs_inf(end)+freq_inf_th),Nd,n_points));
    RE = cat(3,RE,RE_f);
    IM = cat(3,IM,IM_f);
end
NY=(RE(:,:)+1i*IM(:,:));
%%%%%%%%%%%%%%%%%%%%%
NY_L2=PlotLogNyq(NY,n,freqs_inf,mult);
plot(NY_L2{1,1},XPlot.c,'LineWidth',Lw)
for ii = 1:n_freqs_inf
    plot(NY_L2{1,2*ii},XPlot.c,'LineWidth',Lw,'Linestyle',':')
    plot([real(NY_L2{1,2*ii}(end)) real(NY_L2{1,2*ii+1}(1))],[imag(NY_L2{1,2*ii}(end)) imag(NY_L2{1,2*ii+1}(1))],'Color',XPlot.c,'Linewidth',XPlot.Lw)
    plot(NY_L2{1,2*ii+1},XPlot.c,'LineWidth',Lw)
end
% plot for w<0 and complete
compr_fact = .98;
NY_1=NY_L2{1,1};
circ_inf=CompleteNyqPlot(NY_1(1),XPlot,compr_fact);
if XPlot.NgPlot
    plot(compr_fact*conj(NY_L2{1,1}),XPlot.c,'LineWidth',Lw,'Linestyle','--')
    for ii = 1:n_freqs_inf
        plot(compr_fact*conj(NY_L2{1,2*ii}),XPlot.c,'LineWidth',Lw,'Linestyle',':')
        plot(compr_fact*[real(NY_L2{1,2*ii}(end)) real(NY_L2{1,2*ii+1}(1))],-compr_fact*[imag(NY_L2{1,2*ii}(end)) imag(NY_L2{1,2*ii+1}(1))],...
            'Color',XPlot.c,'Linewidth',XPlot.Lw,'Linestyle','--')
        plot(compr_fact*conj(NY_L2{1,2*ii+1}),XPlot.c,'LineWidth',Lw,'Linestyle','--')
    end
end
% Arrows
pp=[150 480 600 700 850];
NY_L2=transpose(cell2mat(NY_L2));
for jj=pp % arrows for w>0
    plotArrow(real(NY_L2(jj)),imag(NY_L2(jj)),real(NY_L2(jj+1)),imag(NY_L2(jj+1)),0.1,0.1,'r')
end
if XPlot.NgPlot
    pp=[150 480 600 700];
    for jj=pp % arrows for w<0
        plotArrow(compr_fact*real(NY_L2(jj)),-compr_fact*imag(NY_L2(jj)),compr_fact*real(NY_L2(jj-1)),-compr_fact*imag(NY_L2(jj-1)),0.1,0.1,'r')
    end
end
grid on
title('Closed Logarithmic Nyquist Plot')
xlabel('Real aXPlotis')
ylabel('Imaginary aXPlotis')

openloop_rhp_poles = sum(real(sys_poles)>0);
fprintf(1,'Number of open-loop poles in RHP: %i\n',openloop_rhp_poles);
NY_L2=[NY_L2;transpose(circ_inf(round(end/2):end))];
[ncirc,npoles_on_im_aXPlotis] = CountEncirc(conj(NY_L2),NY_L2);
if npoles_on_im_aXPlotis > 0
    fprintf(1,'%d closed-loop pole(s) on the im-aXPlotis.\n',npoles_on_im_aXPlotis);
    fprintf(1,'Graph goes through the -1 point,\n');
    fprintf(1,'encirclement counting infeasible.\n');
else
    fprintf(1,'Number of net encirclements around the -1 point:   %i\n',ncirc);
    closedloop_rhp_poles=ncirc+openloop_rhp_poles;
    fprintf(1,'=> Number of closed-loop poles in RHP:   %i\n',...
        closedloop_rhp_poles);
    if closedloop_rhp_poles > 0
        fprintf(1,'=> Closed-loop system is unstable\n');
    else
        fprintf(1,'and no closed-loop poles on Im-aXPlotis\n=> Closed-loop system is asymptotically stable\n');
    end
end

%**************************************************************************
function h=CompleteNyqPlot(NY1,XPlot,compr_fact)
Np_origin=length(find(pole(XPlot.sys)==0)); % Number of poles at the origin
vv=(angle(conj(NY1)):-0.02:angle(conj(NY1))-Np_origin*pi);
delta_spiral=0.003*(Np_origin>2);
xx=max(2,abs(NY1))*(1-delta_spiral*vv).*exp(1i*vv);
if XPlot.NgPlot
    plot(xx,'Color',XPlot.c,'Linewidth',XPlot.Lw,'Linestyle',':')
    plot([real(xx(end)) real(NY1)],[imag(xx(end)) imag(NY1)],'Color',XPlot.c,'Linewidth',XPlot.Lw,'Linestyle','-')
    plot([real(xx(1)) compr_fact*real(NY1)],[imag(xx(1)) -compr_fact*imag(NY1)],'Color',XPlot.c,'Linewidth',XPlot.Lw,'Linestyle','--')
    Nr_arrow=2*Np_origin; % Number of arrows
    if Nr_arrow>0
        for jj=[1:Nr_arrow]
            pp=jj*round(size(xx,2)/(Nr_arrow+1));
            plotArrow(real(xx(pp)),imag(xx(pp)),real(xx(pp+1)),imag(xx(pp+1)),0.1,0.1,'r')
        end
    end
end
h=xx;
return

%**************************************************************************
function circ=ImaginaryPoleCircle(init_point,end_point_abs,mult)
vv=(angle(init_point):-0.02:angle(init_point)-mult*pi);
delta_spiral=0.003*(mult>2);
circ=linspace(abs(init_point),end_point_abs,length(vv)).*(1-delta_spiral*vv).*exp(1i*vv);
return

%**************************************************************************
function plotArrow(XPlot0,y0,XPlot,y,rXPlot,ry,c)
% Plots an arrow at (XPlot,y) along vector dXPlot=XPlot-XPlot0, dy=y-y0. The length of the arrow
% is rXPlot and ry in the directions XPlot and y. The arrow color is c.
Lw=1.2;         % Linewidth
dXPlot=(XPlot-XPlot0)/rXPlot;
dy=(y-y0)/ry;
dXPloty=sqrt(dXPlot^2+dy^2);
fXPlot=-dXPlot/dXPloty;
fy=-dy/dXPloty;
rotpiu=[rXPlot,0;0,ry]*[cos(pi/6), -sin(pi/6); sin(pi/6), cos(pi/6) ]*[fXPlot,fy]';
rotmeno=[rXPlot,0;0,ry]*[cos(pi/6), sin(pi/6); -sin(pi/6), cos(pi/6) ]*[fXPlot,fy]';
plot([XPlot0,XPlot],[y0,y],'-')
plot([XPlot,XPlot+rotpiu(1)],[y,y+rotpiu(2)],c,'LineWidth',Lw)
plot([XPlot,XPlot+rotmeno(1)],[y,y+rotmeno(2)],c,'LineWidth',Lw)
return

%**************************************************************************
function Ln=PlotLogNyq(NY,n,freqs_inf,mult)
RE=real(NY);
IM=imag(NY);
MD=sqrt(RE.^2+IM.^2);
PH=atan2(IM,RE);
MD2=MD.^(log10(n));

ind=find(MD2>1);
MD2(ind)=2-1./MD2(ind);

n_freqs_inf = length(freqs_inf);
n_points = 1000/(n_freqs_inf+1);
Ln = cell(1,2*n_freqs_inf+1);

for ii = 0:n_freqs_inf-1
    Ln{1,1+2*ii}=MD2(1+ii*n_points:(ii+1)*n_points).*exp(1i*PH(1+ii*n_points:(ii+1)*n_points));
    circ_i=ImaginaryPoleCircle(Ln{1,1+2*ii}(end),MD2(1+(ii+1)*n_points),mult(ii+1));
    Ln{1,2+2*ii}=circ_i;
end
Ln{1,end}=MD2(1+n_freqs_inf*n_points:end).*exp(1i*PH(1+n_freqs_inf*n_points:end));
return

%**************************************************************************
function [ncirc,npoles_on_im_aXPlotis] = CountEncirc(zmirr,zmain)

% Counts net encirclements around the point -1.
% An encirclement is counted as positive if the direction
% is clockwise. This follows Belanger (1995):
% "Control Engeering", Saunders College Publishing,
% pp 206 - 208.
%
% Bugs fix and improvements made in Feb. 09:
% The function now also counts poles on the im-axis for the
% closed-loop system, if any. If such poles exist, this
% corresponds to the graph going through -1.
% Encirclement counting is then impossible and is
% disabled.
% Another (small) bug fiXPloted and improvements made in June 09

eps=1e-6;
zmirr(1:end)=zmirr(end:-1:1);
zmirr=zmirr(2:end-1);
zall=[zmirr;zmain;zmirr(1)];
if abs(imag(zall(1))) < eps
    zall=[zall;zall(2)];
end
ncirc=0;
npoles_on_im_aXPlotis=0;
z3=zall(end);
for k=3:length(zall)
    z4=z3;
    z1=zall(k);z2=zall(k-1); z3=zall(k-2);
    abz1=abs(z1+1);abz2=abs(z2+1); abz3=abs(z3+1);
    zre1=real(z1); zre2=real(z2);
    %   Checking if graph is too close to -1:
    dl1= fromline2minusone(z1,z2);
    dl2= fromline2minusone(z2,z3);
    dl3= fromline2minusone(z3,z4);
    closest_now = abz1 > abz2 && abz3 > abz2;
    if closest_now && min([dl1 dl2 dl3]) < 1e-5 ...
            && min([abz1 abz2 abz3])< 0.001
        npoles_on_im_aXPlotis=npoles_on_im_aXPlotis+1;
    end
    %   end checking if graph is too close to -1.
    
    % Only checking for Re aXPlotis crossings to the left
    % of minus 0.9 to avoid unnecessary work:
    if zre1 < -0.9
        zim1=imag(z1); zim2=imag(z2); zim3=imag(z3);
        if zim1*zim2 < 0
            % Interpolation to find real z value at crossing:
            delta12=(real(z1)-real(z2))*abs(imag(z2))/(abs(imag(z1))+abs(imag(z2)));
            realcross=real(z2)+delta12;
        end
        if zim1 > eps && zim2 < -eps
            if realcross < -1
                ncirc=ncirc+1;% ncirc,z1, z2, z3
            end
        elseif zim1 < -eps && zim2 > eps
            if realcross < -1
                ncirc=ncirc-1;% ncirc,z1, z2, z3
            end
        elseif abs(zim2) < eps && zim1 > 0 && zim3 < 0
            if real(z2) < -1
                ncirc=ncirc+1;% ncirc,z1, z2, z3
            end
        elseif abs(zim2) < eps && zim3 > 0 && zim1 < 0
            if real(z2) < -1
                ncirc=ncirc-1;% % ncirc,z1, z2, z3
            end
        else
        end
    end % real(z1) < -0.99
end
return

% ********************************
function [distline_to_minus1] = fromline2minusone(z1,z2)
% Calcucates the min. distance from the point -1 to the line z1,z2
v=z2-z1;
v=imag(v)-1i*real(v);
v=v/abs(z2-z1);
r=z1+1;
d=dot([real(v) imag(v)],[real(r) imag(r)]');
distline_to_minus1=abs(d);
return