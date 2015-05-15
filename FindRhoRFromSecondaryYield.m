function output = FindRhoRFromSecondaryYield(...
    Y1,Y1unc, Y2n, Y2nunc, Y2p, Y2punc, Ti, Te, rho, Zf, Af, fi, HSorUni, varargin)
% Find rhoR that satisfies conditions set by Y2p|Y2n, Y1, Te, Ti, rho
% Inputs:  Y1 = primary yield, Nx1
%          Y1unc = uncertainty in primary yield, Nx1
%          Y2n = secondary neutron yield (0 if unused), Nx1
%          Y2nunc = secondary neutron yield uncertainty, Nx1
%          Y2p = secondary proton yield (0 if unused), Nx1
%          Y2punc = secondary proton yield uncertainty, Nx1
%          Ti = Ion temperature (keV), Nx1
%          Te = electron temperature (keV), Nx1
%          rho = plasma density at BT (g/cc), Nx1
%          Zf = charge state of plasma ions, NxM
%          Af = Mass of plasma ions (amu), NxM
%          fi = number fraction of plasma ions, NxM
%          HSorUni = flag, return hotspot (0) or uniform model (1)
% Outputs: output.rhoR_Y2n = fuel rhoR (mg/cc) from DT-n
%          output.rhoR_Y2p = fuel rhoR (mg/cc) from D3He-p
%          output.rhoR_Y2n_range = [upper, lower] uncertainty range (mg/cc) from DT-n
%          output.rhoR_Y2p_range = [upper, lower] uncertainty range (mg/cc) from D3He-p

%% Test Values
%{
Y1 = 1;
Y1unc = 0;

Y2n = 0;
Y2nunc = 0.1*Y2n;
Y2n_use = 0;

Y2p = 4.95e-4;
Y2punc = 7.68e-5;
Y2p_use = 1;

Ti = 2.6;
Te = 1.1;

rho = 0.21; %g/cm^3
rhounc = 0.09;

Zf = repmat([1,2],[length(Ti),1]);
Af = repmat([2,3],[length(Ti),1]);
fi = repmat([1,0],[length(Ti),1]);
%}

%% Stopping Power parameters
TeCorrection = 1;
quiet = 0; %suppress messages
makeFigures = 1; % make figures

%% Deal with variable arguments
if ~isempty(varargin),
    if mod(length(varargin),2)~=0,
        error('Extra arguments must be declared in pairs ''VariableName'', value');
    end
    
    for i = 1:length(varargin)/2,
        if isstr(varargin{i*2-1});
            arg = varargin{i*2-1};
            val = varargin{i*2};
            if isstr(val),
                eval([arg,'=',val,';'])
            else
                eval([arg,'=',num2str(val),';'])
            end
        end
    end
end

%% Constants (cgs)
e = 4.80320425*10^-10;  %statcoulombs, (erg*cm)^(1/2)
c = 2.99792458*10^10;   %cm/sec
me = 9.10938291*10^-28; %g
mp = 1.67262178*10^-24; %g
mpc2 = 938.272046;      %MeV, proton rest mass
kB = 1.3806488*10^-16;  %erg/K
kBkeV = 1.60217646*10^-9;%erg/keV
Na = 6.02214129*10^23;  %mol^-1
hbar = 1.05457148*10^-27; %erg*sec

%% Correctly shape inputs; Generate Plasma
N = length(Ti);
Ti = reshape(Ti,[N,1]);
Te = reshape(Te,[N,1]);
rho = reshape(rho,[N,1]);

M = size(Zf,2);
Zf = reshape(Zf,[N,M]);
Af = reshape(Af,[N,M]);
fi = reshape(fi,[N,M]);

fi = fi./repmat(sum(fi,2),[1,M]);
f_mass = fi.*Af./repmat(sum(fi.*Af,2),[1,M]);

Dind = ((Zf==1).*(Af==2)==1);
nitot = rho./sum(Af.*fi.*mp,2); %cm^-3
nD = nitot.*fi(Dind);

% Yields
Y1 = reshape(Y1,[N,1]);
Y1unc = reshape(Y1unc,[N,1]);
if Y2n>0, 
    Y2n_use = 1;
    Y2n = reshape(Y2n,[N,1]); 
    Y2nunc = reshape(Y2nunc,[N,1]); 
else
    Y2n_use = 0;
end

if Y2p>0, 
    Y2p_use = 1;
    Y2p = reshape(Y2p,[N,1]); 
    Y2punc = reshape(Y2punc,[N,1]); 
else
    Y2p_use = 0;
end

%% Calculate the secondary yield curves for these plasmas
sec = FindSecondaryYield(Ti, Te, rho, Zf, Af, fi, HSorUni,...
    'TeCorrection',TeCorrection,'quiet', quiet);

PtotDTn = sec.Y2n;
PtotD3Hep = sec.Y2p;
rhoRt = sec.rhoR_Y2n;
rhoR3He = sec.rhoR_Y2p;

% output.test1 = PtotDTn;
% output.test2 = PtotD3Hep;
% output.test3 = rhoRt;
% output.test4 = rhoR3He;

%% Calculate fuel rhoR for each species
if Y2n_use,
    fuelRhoR_n_index = sum(PtotDTn<=repmat((Y2n./Y1),[1,size(PtotDTn,2)]),2);
    fuelRhoR_n_unc_index{1} = sum(PtotDTn<=repmat((Y2n./Y1).*(1+sqrt((Y2nunc./Y2n).^2+(Y1unc./Y1).^2)),[1,size(PtotDTn,2)]),2);
    fuelRhoR_n_unc_index{2} = sum(PtotDTn<=repmat((Y2n./Y1).*(1-sqrt((Y2nunc./Y2n).^2+(Y1unc./Y1).^2)),[1,size(PtotDTn,2)]),2);
    fuelRhoR_n = diag(rhoRt(:,fuelRhoR_n_index));
    fuelRhoR_n_unc(:,1) = diag(rhoRt(:,fuelRhoR_n_unc_index{1}));
    fuelRhoR_n_unc(:,2) = diag(rhoRt(:,fuelRhoR_n_unc_index{2}));

    if ~quiet,
        fprintf('\nfrom Secondary DT-neutrons:\n');
        for i = 1:length(fuelRhoR_n),
            if length(fuelRhoR_n)>1, fprintf('Plasma #%g:\n',i); end
            fprintf('Total rhoR in fuel (mg/cm2): %g\n',fuelRhoR_n(i)*1000);
            fprintf('       uncertainty: [+%g,%g]\n',(fuelRhoR_n_unc(i,:)-fuelRhoR_n(i))*1000);
            fprintf('Deuterium rhoR (mg/cm2): %g\n',fuelRhoR_n(i)*1000.*f_mass(i,Dind(i,:)));
        end
        fprintf('(note:  uncertainty due to uncertainty in yield only!)\n');
    end
else
    fuelRhoR_n = zeros([N,1]);
    fuelRhoR_n_unc = zeros([N,2]);
end

if Y2p_use,
    fuelRhoR_p_index = sum(PtotD3Hep<=repmat((Y2p./Y1),[1,size(PtotD3Hep,2)]),2);
    fuelRhoR_p_unc_index{1} = sum(PtotD3Hep<=repmat((Y2p./Y1).*(1+sqrt((Y2punc./Y2p).^2+(Y1unc./Y1).^2)),[1,size(PtotD3Hep,2)]),2);
    fuelRhoR_p_unc_index{2} = sum(PtotD3Hep<=repmat((Y2p./Y1).*(1-sqrt((Y2punc./Y2p).^2+(Y1unc./Y1).^2)),[1,size(PtotD3Hep,2)]),2);
    fuelRhoR_p = diag(rhoR3He(:,fuelRhoR_p_index));
    fuelRhoR_p_unc(:,1) = diag(rhoR3He(:,fuelRhoR_p_unc_index{1}));
    fuelRhoR_p_unc(:,2) = diag(rhoR3He(:,fuelRhoR_p_unc_index{2}));

    if ~quiet,
        fprintf('\nfrom Secondary D3He-protons:\n');
        for i = 1:length(fuelRhoR_p),
            if length(fuelRhoR_p)>1, fprintf('Plasma #%g:\n',i); end
            fprintf('Total rhoR in fuel (mg/cm2): %g\n',fuelRhoR_p(i)*1000);
            fprintf('       uncertainty: [+%g,%g]\n',(fuelRhoR_p_unc(i,:)-fuelRhoR_p(i))*1000);
            fprintf('Deuterium rhoR (mg/cm2): %g\n',fuelRhoR_p(i)*1000.*f_mass(i,Dind(i,:)));
        end
        fprintf('(note:  uncertainty due to uncertainty in yield only!)\n');
    end
else
    fuelRhoR_p = zeros([N,1]);
    fuelRhoR_p_unc = zeros([N,2]);
end
fprintf('\n');

output.rhoR_Y2n = fuelRhoR_n*1000;
output.rhoR_Y2p = fuelRhoR_p*1000;
output.rhoR_Y2n_range = [(fuelRhoR_n-fuelRhoR_n_unc(:,2))*1000, (fuelRhoR_n_unc(:,1)-fuelRhoR_n)*1000];
output.rhoR_Y2p_range = [(fuelRhoR_p-fuelRhoR_p_unc(:,2))*1000, (fuelRhoR_p_unc(:,1)-fuelRhoR_p)*1000];

%% Make plots of the ranges & measurements
if makeFigures,
    figure; 
    colors = {'r','b','g','m','k','c','y'};
    for i = 1:N,
        if Y2n_use, 
            loglog(rhoRt(i,:),PtotDTn(i,:),colors{i},'LineWidth',2); 
            hold on;
            plot(fuelRhoR_n_unc,[1,1]*Y2n(i)/Y1(i),[colors{i},':'],'LineWidth',1);
            plot([1,1]*fuelRhoR_n,Y2n(i)/Y1(i)*[1+sqrt((Y2nunc(i)/Y2n(i))^2+(Y1unc(i)/Y1(i))^2),1-sqrt((Y2nunc(i)/Y2n(i))^2+(Y1unc(i)/Y1(i))^2)],...
                [colors{i},':'],'LineWidth',1);
            plot(fuelRhoR_n(i),Y2n(i)/Y1(i),[colors{i},'x'],'LineWidth',2,'MarkerSize',12);
        end
        if Y2p_use, 
            loglog(rhoR3He(i,:),PtotD3Hep(i,:),colors{i},'LineWidth',2); 
            hold on;
            plot(fuelRhoR_p_unc,[1,1]*Y2p(i)/Y1(i),[colors{i},':'],'LineWidth',1);
            plot([1,1]*fuelRhoR_p,Y2p(i)/Y1(i)*[1+sqrt((Y2punc(i)/Y2p(i))^2+(Y1unc(i)/Y1(i))^2),1-sqrt((Y2punc(i)/Y2p(i))^2+(Y1unc(i)/Y1(i))^2)],...
                [colors{i},':'],'LineWidth',1);
            plot(fuelRhoR_p(i),Y2p(i)/Y1(i),[colors{i},'x'],'LineWidth',2,'MarkerSize',12);
        end
        set(gca,'FontSize',14);
        if HSorUni, title('Hotspot Model');
        else title('Uniform Model'); end
        
        xlabel('fuel RhoR (g/cc)');
        ylabel('Y2/Y1');
        box;
    end
end
