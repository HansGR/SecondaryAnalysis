function output = FindSecondaryYield(Ti, Te, rho, Zf, Af, fi, HSorUni, varargin)
% Find the Y2/Y1 curves for a plasma with Te, Ti, rho, composition
% Inputs:  Ti = Ion temperature (keV), Nx1
%          Te = electron temperature (keV), Nx1
%          rho = plasma density at BT (g/cc), Nx1
%          Zf = charge state of plasma ions, NxM
%          Af = Mass of plasma ions (amu), NxM
%          fi = number fraction of plasma ions, NxM
%          HSorUni = flag, return hotspot (0) or uniform model (1)
% Outputs: for each input set,
%          output.rhoR_Y2n = vector rhoR axis for Y2n 
%          output.Y2n      = vector Y2n/Y1
%          output.rhoR_Y2p = vector rhoR axis for Y2p
%          output.Y2p      = vector Y2p/Y1

%% Test Values
%{
Ti = 1:1:10;
Te = Ti;

rho = 0.21*ones(size(Ti)); %g/cm^3

Zf = repmat([1,2],[length(Ti),1]);
Af = repmat([2,3],[length(Ti),1]);
fi = repmat([1,0],[length(Ti),1]);

HSorUni = 0;
%}

%% Stopping Power parameters
TeCorrection = 1;
mimicFredrick = 0;
%quiet = 1; %suppress messages
electronsOnly = 0;  % only use electron stopping power in the calculation
output_escapeEnergy = 0;  % include the (mean) escaping energy of the tritons/3He in the output
quiet = 1; 
parallel = 0;

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

fi = fi./repmat(sum(fi,2),[1,size(fi,2)]);
f_mass = fi.*Af./repmat(sum(fi.*Af,2),[1,size(fi,2)]);

Dind = ((Zf==1).*(Af==2)==1);
if Dind==0,
    error('Fuel must contain deuterium!')
end

nitot = rho./sum(Af.*fi.*mp,2); %cm^-3
nD = nitot.*fi(Dind);

%% Establish x axes
%path(path,'E:\Data\Code\Matlab\LiPetrasso Stopping Power');

E0t = 1.01;  %MeV triton
Z0t = 1*ones(N,1);
A0t = 3*ones(N,1);

E0_3He = 0.82; %MeV 3He
Z0_3He = 2*ones(N,1);
A0_3He = 3*ones(N,1);

Et = zeros(N,1);
Et(:,1) = E0t;

E3He = zeros(N,1);
E3He(:,1) = E0_3He;

% rangeFactor = 3;
% stppwrt = LiPetrassoStoppingPower(Et,Z0t,A0t,Zf,Af,fi,nitot,...
%         Ti,Te,'MeV/um',TeCorrection,1,mimicFredrick);
% if ~electronsOnly,
%     stppwrt = sum(stppwrt,2);
% else
%     stppwrt = stppwrt(:,end);
% end
% expectedRangeT = E0t./stppwrt*rangeFactor;
% 
% stppwr3He = LiPetrassoStoppingPower(E3He,Z0_3He,A0_3He,Zf,Af,fi,nitot,...
%         Ti,Te,'MeV/um',TeCorrection,1,mimicFredrick);
% if ~electronsOnly,
%     stppwr3He = sum(stppwr3He,2);
% else
%     stppwr3He = stppwr3He(:,end);
% end
% expectedRange3He = E0_3He./stppwr3He*rangeFactor;
rangeFactor = 1.5;
dE = 0.01;
Etest = dE:dE:1;
Etest = Etest';

expectedRangeT = zeros([N,1]);
expectedRange3He = zeros([N,1]);
for i = 1:N,
    stppwrt_temp = LiPetrassoStoppingPower(E0t*Etest,Z0t(i)*ones(size(Etest)),A0t(i)*ones(size(Etest)),...
        repmat(Zf(i,:),[length(Etest),1]),repmat(Af(i,:),[length(Etest),1]),repmat(fi(i,:),[length(Etest),1]),...
        nitot(i)*ones(size(Etest)),Ti(i)*ones(size(Etest)),Te(i)*ones(size(Etest)),...
        'MeV/um',TeCorrection,1,mimicFredrick);
    if ~electronsOnly,
        stppwrt(i,:) = sum(stppwrt_temp,2);
    else
        stppwrt(i,:) = stppwrt_temp(:,end);
    end
    beg_ind_t(i) = find(stppwrt(i,:)>0,1,'first');
    expectedRangeT(i) = sum(1./stppwrt(i,beg_ind_t(i):end))*E0t*dE*rangeFactor; %microns
    
    stppwr3He_temp = LiPetrassoStoppingPower(E0_3He*Etest,Z0_3He(i)*ones(size(Etest)),A0_3He(i)*ones(size(Etest)),...
        repmat(Zf(i,:),[length(Etest),1]),repmat(Af(i,:),[length(Etest),1]),repmat(fi(i,:),[length(Etest),1]),...
        nitot(i)*ones(size(Etest)),Ti(i)*ones(size(Etest)),Te(i)*ones(size(Etest)),...
        'MeV/um',TeCorrection,0,mimicFredrick);
    if ~electronsOnly,
        stppwr3He(i,:) = sum(stppwr3He_temp,2);
    else
        stppwr3He(i,:) = stppwr3He_temp(:,end);
    end
    beg_ind_3He(i) = find(stppwr3He(i,:)>0,1,'first');
    expectedRange3He(i) = sum(1./stppwr3He(i,beg_ind_3He(i):end))*E0_3He*dE*rangeFactor; %microns
end

M = 2000;                   %number of steps
dxt = expectedRangeT/M;     %um
dx3He = expectedRange3He/M; %um

xt = repmat(1:M,[N,1]);       
xt = xt.*repmat(dxt,[1,M]);    %um
x3He = repmat(1:M,[N,1]);
x3He = x3He.*repmat(dx3He,[1,M]); %um

rhoRt = (xt*1e-4).*repmat(rho,[1,M]); %g/cm^2
rhoR3He = (x3He*1e-4).*repmat(rho,[1,M]); %g/cm^2

rhoRt_d = (xt*1e-4).*repmat(nD.*2.*mp,[1,M]); %g/cm^2, deuterium only
rhoR3He_d = (x3He*1e-4).*repmat(nD.*2.*mp,[1,M]); %g/cm^2, deuterium only

%% Energy; Stopping Power
Et = zeros(N,M+1);
Et(:,1) = E0t;

E3He = zeros(N,M+1);
E3He(:,1) = E0_3He;

if ~quiet, fprintf('Ranging particles...'); end
for i = 1:M,
    %fprintf('%g\n',i);
    logical = Et(:,i)>0;
    if sum(logical)>0,
        stppwrt = LiPetrassoStoppingPower(Et(logical,i),Z0t(logical),A0t(logical),...
            Zf(logical,:),Af(logical,:),fi(logical,:),nitot(logical),...
            Ti(logical),Te(logical),'MeV/um',TeCorrection,1,mimicFredrick);
        if ~electronsOnly,
            stppwrt = sum(stppwrt,2);
        else
            stppwrt = stppwrt(:,end);
        end
        Et(logical,i+1) = Et(logical,i) - stppwrt.*dxt(logical);
        Et(logical,i+1) = Et(logical,i+1).*(stppwrt>0).*(Et(logical,i)*1000>Te(logical));
    end
    
    logical_3He = E3He(:,i)>0;
    if sum(logical_3He)>0,
        stppwr3He = LiPetrassoStoppingPower(E3He(logical_3He,i),Z0_3He(logical_3He),A0_3He(logical_3He),...
            Zf(logical_3He,:),Af(logical_3He,:),fi(logical_3He,:),nitot(logical_3He),...
            Ti(logical_3He),Te(logical_3He),'MeV/um',TeCorrection,1,mimicFredrick);
        if ~electronsOnly,
            stppwr3He = sum(stppwr3He,2);
        else
            stppwr3He = stppwr3He(:,end);
        end
        E3He(logical_3He,i+1) = E3He(logical_3He,i) - stppwr3He.*dx3He(logical_3He);
        E3He(logical_3He,i+1) = E3He(logical_3He,i+1).*(stppwr3He>0).*(E3He(logical_3He,i)*1000>Te(logical_3He));
    end
end
if ~quiet, fprintf('done\n'); end

Et(Et<0) = 0;
E3He(E3He<0) = 0;

%% Figures
%{
figure;
for i = 11:15, semilogx(xt(i,:),Et(i,:)); hold on; end
title('Triton stopping');

figure; 
for i = 11:15, semilogx(x3He(i,:),E3He(i,:)); hold on; end
title('3He stopping');
%}
%% Cross section
%path(path,'E:\Data\Code\Matlab\HYADES');
Et_cell = (Et(:,1:end-1)+Et(:,2:end))/2;
E3He_cell = (E3He(:,1:end-1)+E3He(:,2:end))/2;

ECOMt = 1000*Et_cell.*repmat((2./(A0t+2)),[1,M]) + repmat(Ti,[1,M]).*repmat((A0t./(A0t+2)),[1,M]); %triton COM energy, keV
ECOM3He = 1000*E3He_cell.*repmat((2./(A0_3He+2)),[1,M]) + repmat(Ti,[1,M]).*repmat((A0_3He./(A0_3He+2)),[1,M]); %3He COM energy, keV

cst = BoschHaleCrossSection(ECOMt,'DT');      %mb = 1e-27 cm^2
cs3He = BoschHaleCrossSection(ECOM3He,'D3He');  %mb = 1e-27 cm^2

%% Figures
%{
figure;
for i = 11:15, plot(xt(i,:),cst(i,:)); hold on; end
title('Triton xs'); xlabel('distance (\mu m)');

figure; 
for i = 11:15, plot(x3He(i,:),cs3He(i,:)); hold on; end
title('3He xs'); xlabel('distance (\mu m)');
%}

%% Probability of fusion:  Hotspot model
ProbDTn = repmat(nD,[1,M]).*cst*1e-27;  %probability of DT-n fusion per cm (dP/dx), plot vs x
ProbDTn(isnan(ProbDTn)) = 0;

ProbD3Hep = repmat(nD,[1,M]).*cs3He*1e-27;  %probability of D3He-p fusion per cm (dP/dx), plot vs x
ProbD3Hep(isnan(ProbD3Hep)) = 0;

%Total probability of 2ndary fusion within a radius (given by rhoRt,rhoR3He)
PtotDTn = cumsum(ProbDTn.*repmat(dxt*1e-4,[1,M]),2);
PtotD3Hep = cumsum(ProbD3Hep.*repmat(dx3He*1e-4,[1,M]),2);

if HSorUni==0,
    output.rhoR_Y2n = rhoRt;
    output.rhoR_Y2p = rhoR3He;
    output.Y2n = PtotDTn;
    output.Y2p = PtotD3Hep;
    
    if output_escapeEnergy,
        output.Et = Et_cell;
        output.E3He = E3He_cell;
    end
    % Calculate the fuel rhoR from Y2n and/or Y2p (HOTSPOT)
%     if ~quiet, fprintf('HOTSPOT MODEL:\n'); end
%     if Y2n_use,
%         fuelRhoR_n_index = sum(PtotDTn<=repmat((Y2n./Y1),[1,size(PtotDTn,2)]),2);
%         fuelRhoR_n_unc_index{1} = sum(PtotDTn<=repmat((Y2n./Y1).*(1+sqrt((Y2nunc./Y2n).^2+(Y1unc./Y1).^2)),[1,size(PtotDTn,2)]),2);
%         fuelRhoR_n_unc_index{2} = sum(PtotDTn<=repmat((Y2n./Y1).*(1-sqrt((Y2nunc./Y2n).^2+(Y1unc./Y1).^2)),[1,size(PtotDTn,2)]),2);
%         fuelRhoR_n = diag(rhoRt(:,fuelRhoR_n_index));
%         fuelRhoR_n_unc(:,1) = diag(rhoRt(:,fuelRhoR_n_unc_index{1}));
%         fuelRhoR_n_unc(:,2) = diag(rhoRt(:,fuelRhoR_n_unc_index{2}));
% 
%         if ~quiet,
%             fprintf('from Secondary DT-neutrons:\n');
%             fprintf('Total rhoR in fuel (mg/cm2): %g\n',fuelRhoR_n*1000);
%             fprintf('       uncertainty: %g\n',abs(diff(fuelRhoR_n_unc)/2)*1000);
%             fprintf('Deuterium rhoR (mg/cm2): %g\n',fuelRhoR_n*1000.*f_mass(Dind));
%             fprintf('(note:  uncertainty due to uncertainty in yield only!)\n');
%         end
%     else
%         fuelRhoR_n = zeros([N,1]);
%         fuelRhoR_n_unc = zeros([N,2]);
%     end
% 
%     if Y2p_use,
%         fuelRhoR_p_index = sum(PtotD3Hep<=repmat((Y2p./Y1),[1,size(PtotD3Hep,2)]),2);
%         fuelRhoR_p_unc_index{1} = sum(PtotD3Hep<=repmat((Y2p./Y1).*(1+sqrt((Y2punc./Y2p).^2+(Y1unc./Y1).^2)),[1,size(PtotD3Hep,2)]),2);
%         fuelRhoR_p_unc_index{2} = sum(PtotD3Hep<=repmat((Y2p./Y1).*(1-sqrt((Y2punc./Y2p).^2+(Y1unc./Y1).^2)),[1,size(PtotD3Hep,2)]),2);
%         fuelRhoR_p = diag(rhoR3He(:,fuelRhoR_p_index));
%         fuelRhoR_p_unc(:,1) = diag(rhoR3He(:,fuelRhoR_p_unc_index{1}));
%         fuelRhoR_p_unc(:,2) = diag(rhoR3He(:,fuelRhoR_p_unc_index{2}));
%         
%         if ~quiet,
%             fprintf('from Secondary D3He-protons:\n');
%             fprintf('Total rhoR in fuel (mg/cm2): %g\n',fuelRhoR_p*1000);
%             fprintf('       uncertainty: %g\n',abs((fuelRhoR_p_unc(:,1)-fuelRhoR_p_unc(:,2))/2)*1000);
%             fprintf('Deuterium rhoR (mg/cm2): %g\n',fuelRhoR_p*1000.*f_mass(Dind));
%             fprintf('(note:  uncertainty due to uncertainty in yield only!)\n');
%         end
%     else
%         fuelRhoR_p = zeros([N,1]);
%         fuelRhoR_p_unc = zeros([N,2]);
%     end
%     fprintf('\n');
% 
%     rhoR = [fuelRhoR_n*1000,fuelRhoR_p*1000];
%     rhoRunc_hi = [(fuelRhoR_n_unc(:,1)-fuelRhoR_n)*1000,...
%         (fuelRhoR_p_unc(:,1)-fuelRhoR_p)*1000];
%     rhoRunc_lo = [(fuelRhoR_n-fuelRhoR_n_unc(:,2))*1000,...
%         (fuelRhoR_p-fuelRhoR_p_unc(:,2))*1000];
% 
%     % Make plots of the ranges & measurements
%     figure; hold on;
%     colors = {'r','b','g','m','k','c'};
%     for i = 1:N,
%         if Y2n_use, 
%             plot(rhoRt(i,:),PtotDTn(i,:),colors{i},'LineWidth',2); 
%             plot(fuelRhoR_n(i),Y2n(i)/Y1(i),[colors{i},'x'],'LineWidth',2,'MarkerSize',12);
%         end
%         if Y2p_use, 
%             plot(rhoR3He(i,:),PtotD3Hep(i,:),colors{i},'LineWidth',2); 
%             plot(fuelRhoR_p(i),Y2p(i)/Y1(i),[colors{i},'x'],'LineWidth',2,'MarkerSize',12);
%         end
%         title('Hotspot Model');
%         xlabel('fuel RhoR (g/cc)','FontSize',14);
%         ylabel('Y2/Y1','FontSize',14);
%         set(gca,'FontSize',14);
%         box;
%     end
    
else
    %% Probability of fusion: Uniform model
    %The uniform model results can be obtained from an appropriate numerical
    %average of the hotspot results; basically, averaging the prob of 2ndary
    %fusion over 4 pi from each radial point, then averaging over the sphere.

    %Extend rhoRt, rhoR3He to include much larger values.  This is OK because
    %the hotspot probabilities should have all saturated already.  If they
    %haven't this won't be correct.
    numExtend = 10;
    maxExtend = 1; %g/cm^2

    delta_t = (maxExtend./rhoRt(:,M)).^(1/numExtend);
    delta_t_mult = repmat(delta_t,[1,numExtend]).^repmat([1:numExtend],[N,1]);
    rhoRt_ext = rhoRt;
    rhoRt_ext(:,M+1:M+numExtend) = repmat(rhoRt(:,M),[1,numExtend]).*delta_t_mult;

    delta_3He = (maxExtend./rhoR3He(:,M)).^(1/numExtend);
    delta_3He_mult = repmat(delta_3He,[1,numExtend]).^repmat([1:numExtend],[N,1]);
    rhoR3He_ext = rhoR3He;
    rhoR3He_ext(:,M+1:M+numExtend) = repmat(rhoR3He(:,M),[1,numExtend]).*delta_3He_mult; %g/cm^3

    PtotDTn_ext = PtotDTn;
    PtotDTn_ext(:,M+1:M+numExtend) = repmat(PtotDTn(:,M),[1,numExtend]);
    PtotD3Hep_ext = PtotD3Hep;
    PtotD3Hep_ext(:,M+1:M+numExtend) = repmat(PtotD3Hep(:,M),[1,numExtend]);
    
    if output_escapeEnergy,
        Et_ext = Et_cell;
        Et_ext(:,M+1:M+numExtend) = repmat(E3He(:,M),[1,numExtend]);
        E3He_ext = E3He_cell;
        E3He_ext(:,M+1:M+numExtend) = repmat(E3He(:,M),[1,numExtend]);
        
        Et_u = zeros(size(Et_ext));
        E3He_u = zeros(size(E3He_ext));
    end
    
    %define radius, theta points for average
    numR = 100; %number of discrete radii
    numTh = 20; %number of angles

    %matrix dimensions will be: [rhoR_total (M+1+numExtend), r steps (numR), theta (numTh)]
    r = reshape([0:1/(numR-1):1],[1,numR]);
    rmat = repmat(r,[M+numExtend,1,numTh]);

    th = reshape([0:1/(numTh-1):1]*pi,[1,1,numTh]);
    thmat = repmat(th,[M+numExtend,numR,1]);

    rth_factor = -rmat.*cos(thmat)+sqrt(1-rmat.^2.*sin(thmat).^2);

    %calculate Uniform Model probability
    if ~quiet,  fprintf('Calculating Uniform model:\n'); end
    PtotDTn_u = zeros(size(PtotDTn_ext));
    PtotD3Hep_u = zeros(size(PtotD3Hep_ext));

    if parallel,
        if ~output_escapeEnergy, 
            Et_ext = zeros([N,1]); 
            E3He_ext = zeros([N,1]); 
            Et_rth = 0;
            Et_u = zeros([N,1]); 
            E3He_rth = 0;
            E3He_u = zeros([N,1]); 
        end
        
        if isempty(gcp('nocreate'))
            parpool;
        end

        parfor_progress(N);
        parfor i = 1:N,
            rhoRt_rth = repmat(reshape(rhoRt_ext(i,:),[M+numExtend,1]),[1,numR,numTh]).*rth_factor;
            PtotDTn_rth = interp1([0,rhoRt_ext(i,:)],[0,PtotDTn_ext(i,:)],rhoRt_rth,'pchip',max(PtotDTn(i,:)));
            PtotDTn_u(i,:) = sum(sum(PtotDTn_rth.*rmat.^2.*sin(thmat),3),2)*(pi/(numTh-1))*(1/(numR-1))...
                /(sum(sum(rmat(1,:,:).^2.*sin(thmat(1,:,:)),3),2)*(pi/(numTh-1))*(1/(numR-1))); %normalization factor: 2/3, in limit

            rhoR3He_rth = repmat(reshape(rhoR3He_ext(i,:),[M+numExtend,1]),[1,numR,numTh]).*rth_factor;
            PtotD3Hep_rth = interp1([0,rhoR3He_ext(i,:)],[0,PtotD3Hep_ext(i,:)],rhoR3He_rth,'pchip',max(PtotD3Hep(i,:)));
            PtotD3Hep_u(i,:) = sum(sum(PtotD3Hep_rth.*rmat.^2.*sin(thmat),3),2)*(pi/(numTh-1))*(1/(numR-1))...
                /(sum(sum(rmat(1,:,:).^2.*sin(thmat(1,:,:)),3),2)*(pi/(numTh-1))*(1/(numR-1))); %normalization factor: 2/3, in limit

            if output_escapeEnergy,
                % Calculate the mean escaping energy from a plasma with rhoRt
                Et_rth = interp1([0,rhoRt_ext(i,:)],[E0t,Et_ext(i,:)],rhoRt_rth,'pchip',min(Et(i,:)));
                Et_u(i,:) = sum(sum(Et_rth.*rmat.^2.*sin(thmat),3),2)*(pi/(numTh-1))*(1/(numR-1))...
                    /(sum(sum(rmat(1,:,:).^2.*sin(thmat(1,:,:)),3),2)*(pi/(numTh-1))*(1/(numR-1))); %normalization factor: 2/3, in limit

                E3He_rth = interp1([0,rhoRt_ext(i,:)],[E0_3He,E3He_ext(i,:)],rhoRt_rth,'pchip',min(E3He(i,:)));
                E3He_u(i,:) = sum(sum(E3He_rth.*rmat.^2.*sin(thmat),3),2)*(pi/(numTh-1))*(1/(numR-1))...
                    /(sum(sum(rmat(1,:,:).^2.*sin(thmat(1,:,:)),3),2)*(pi/(numTh-1))*(1/(numR-1))); %normalization factor: 2/3, in limit
            end
            parfor_progress;
        end
        parfor_progress(0);
    else
        for i = 1:N,
            if ~quiet, fprintf('%g of %g\n',i,N);    end
            rhoRt_rth = repmat(reshape(rhoRt_ext(i,:),[M+numExtend,1]),[1,numR,numTh]).*rth_factor;
            PtotDTn_rth = interp1([0,rhoRt_ext(i,:)],[0,PtotDTn_ext(i,:)],rhoRt_rth,'pchip',max(PtotDTn(i,:)));
            PtotDTn_u(i,:) = sum(sum(PtotDTn_rth.*rmat.^2.*sin(thmat),3),2)*(pi/(numTh-1))*(1/(numR-1))...
                /(sum(sum(rmat(1,:,:).^2.*sin(thmat(1,:,:)),3),2)*(pi/(numTh-1))*(1/(numR-1))); %normalization factor: 2/3, in limit

            rhoR3He_rth = repmat(reshape(rhoR3He_ext(i,:),[M+numExtend,1]),[1,numR,numTh]).*rth_factor;
            PtotD3Hep_rth = interp1([0,rhoR3He_ext(i,:)],[0,PtotD3Hep_ext(i,:)],rhoR3He_rth,'pchip',max(PtotD3Hep(i,:)));
            PtotD3Hep_u(i,:) = sum(sum(PtotD3Hep_rth.*rmat.^2.*sin(thmat),3),2)*(pi/(numTh-1))*(1/(numR-1))...
                /(sum(sum(rmat(1,:,:).^2.*sin(thmat(1,:,:)),3),2)*(pi/(numTh-1))*(1/(numR-1))); %normalization factor: 2/3, in limit

            if output_escapeEnergy,
                % Calculate the mean escaping energy from a plasma with rhoRt
                Et_rth = interp1([0,rhoRt_ext(i,:)],[E0t,Et_ext(i,:)],rhoRt_rth,'pchip',min(Et(i,:)));
                Et_u(i,:) = sum(sum(Et_rth.*rmat.^2.*sin(thmat),3),2)*(pi/(numTh-1))*(1/(numR-1))...
                    /(sum(sum(rmat(1,:,:).^2.*sin(thmat(1,:,:)),3),2)*(pi/(numTh-1))*(1/(numR-1))); %normalization factor: 2/3, in limit

                E3He_rth = interp1([0,rhoRt_ext(i,:)],[E0_3He,E3He_ext(i,:)],rhoRt_rth,'pchip',min(E3He(i,:)));
                E3He_u(i,:) = sum(sum(E3He_rth.*rmat.^2.*sin(thmat),3),2)*(pi/(numTh-1))*(1/(numR-1))...
                    /(sum(sum(rmat(1,:,:).^2.*sin(thmat(1,:,:)),3),2)*(pi/(numTh-1))*(1/(numR-1))); %normalization factor: 2/3, in limit
            end
        end
    end
    output.rhoR_Y2n = rhoRt_ext;
    output.rhoR_Y2p = rhoR3He_ext;
    output.Y2n = PtotDTn_u;
    output.Y2p = PtotD3Hep_u;
    
    if output_escapeEnergy,
        output.Et = Et_u;
        output.E3He = E3He_u;
    end
    % Calculate the fuel rhoR from Y2n and/or Y2p (UNIFORM)
%     if ~quiet, fprintf('UNIFORM MODEL:\n'); end
%     if Y2n_use,
%         fuelRhoR_n_index = sum(PtotDTn_u<=repmat((Y2n./Y1),[1,size(PtotDTn_u,2)]),2);
%         fuelRhoR_n_unc_index{1} = sum(PtotDTn_u<=repmat((Y2n./Y1).*(1+sqrt((Y2nunc./Y2n).^2+(Y1unc./Y1).^2)),[1,size(PtotDTn_u,2)]),2);
%         fuelRhoR_n_unc_index{2} = sum(PtotDTn_u<=repmat((Y2n./Y1).*(1-sqrt((Y2nunc./Y2n).^2+(Y1unc./Y1).^2)),[1,size(PtotDTn_u,2)]),2);
%         fuelRhoR_n = diag(rhoRt_ext(:,fuelRhoR_n_index));
%         fuelRhoR_n_unc(:,1) = diag(rhoRt_ext(:,fuelRhoR_n_unc_index{1}));
%         fuelRhoR_n_unc(:,2) = diag(rhoRt_ext(:,fuelRhoR_n_unc_index{2}));
%         
%         if ~quiet,
%             fprintf('from Secondary DT-neutrons:\n');
%             fprintf('Total rhoR in fuel (mg/cm2): %g\n',fuelRhoR_n*1000);
%             fprintf('       uncertainty: %g\n',abs(diff(fuelRhoR_n_unc)/2)*1000);
%             fprintf('Deuterium rhoR (mg/cm2): %g\n',fuelRhoR_n*1000.*f_mass(Dind));
%             fprintf('(note:  uncertainty due to uncertainty in yield only!)\n');
%         end
%     else
%         fuelRhoR_n = zeros([N,1]);
%         fuelRhoR_n_unc = zeros([N,2]);
%     end
% 
%     if Y2p_use,
%         fuelRhoR_p_index = sum(PtotD3Hep_u<=repmat((Y2p./Y1),[1,size(PtotD3Hep_u,2)]),2);
%         fuelRhoR_p_unc_index{1} = sum(PtotD3Hep_u<=repmat((Y2p./Y1).*(1+sqrt((Y2punc./Y2p).^2+(Y1unc./Y1).^2)),[1,size(PtotD3Hep_u,2)]),2);
%         fuelRhoR_p_unc_index{2} = sum(PtotD3Hep_u<=repmat((Y2p./Y1).*(1-sqrt((Y2punc./Y2p).^2+(Y1unc./Y1).^2)),[1,size(PtotD3Hep_u,2)]),2);
%         fuelRhoR_p = diag(rhoR3He_ext(:,fuelRhoR_p_index));
%         fuelRhoR_p_unc(:,1) = diag(rhoR3He_ext(:,fuelRhoR_p_unc_index{1}));
%         fuelRhoR_p_unc(:,2) = diag(rhoR3He_ext(:,fuelRhoR_p_unc_index{2}));
% 
%         if ~quiet,
%             fprintf('from Secondary D3He-protons:\n');
%             fprintf('Total rhoR in fuel (mg/cm2): %g\n',fuelRhoR_p*1000);
%             fprintf('       uncertainty: %g\n',abs(diff(fuelRhoR_p_unc)/2)*1000);
%             fprintf('Deuterium rhoR (mg/cm2): %g\n',fuelRhoR_p*1000.*f_mass(Dind));
%             fprintf('(note:  uncertainty due to uncertainty in yield only!)\n');
%         end
%     else
%         fuelRhoR_p = zeros([N,1]);
%         fuelRhoR_p_unc = zeros([N,2]);
%     end
%     fprintf('\n');
% 
%     rhoR = [fuelRhoR_n*1000,fuelRhoR_p*1000];
%     rhoRunc_hi = [(fuelRhoR_n_unc(:,1)-fuelRhoR_n)*1000,...
%         (fuelRhoR_p_unc(:,1)-fuelRhoR_p)*1000];
%     rhoRunc_lo = [(fuelRhoR_n-fuelRhoR_n_unc(:,2))*1000,...
%         (fuelRhoR_p-fuelRhoR_p_unc(:,2))*1000];
%     
%     % Make plots of the ranges & measurements
%     figure; hold on;
%     colors = {'r','b','g','m','k','c','y'};
%     for i = 1:N,
%         if Y2n_use, 
%             plot(rhoRt_ext(i,:),PtotDTn_u(i,:),colors{i},'LineWidth',2); 
%             plot(fuelRhoR_n(i),Y2n(i)/Y1(i),[colors{i},'x'],'LineWidth',2,'MarkerSize',12);
%         end
%         if Y2p_use, 
%             plot(rhoR3He_ext(i,:),PtotD3Hep_u(i,:),colors{i},'LineWidth',2); 
%             plot(fuelRhoR_p(i),Y2p(i)/Y1(i),[colors{i},'x'],'LineWidth',2,'MarkerSize',12);
%         end
%         title('Uniform Model');
%         xlabel('fuel RhoR (g/cc)','FontSize',14);
%         ylabel('Y2/Y1','FontSize',14);
%         set(gca,'FontSize',14);
%         box;
%     end
end