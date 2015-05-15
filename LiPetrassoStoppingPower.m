function output = LiPetrassoStoppingPower(Et,Zt,At,Zf,Af,fi,nitot,Ti,Te,varargin)
%Calculates the stopping power for a test particle in a plasma, using the 
%Li-Petrasso model (PRL 1993).  Accepts a list of entries, length L. Plasma
%contains N species, and has defined for each list entry: fi (fraction of 
%each species), nitot, Ti, and Te. 
%
%Inputs: All vectors, length L unless noted
%  Test ion:
%   Et = test particle energy, MeV
%   Zt = test particle charge 
%   At = test particle atomic number
%  Field plasma properties:
%   fi = field ion fractional number density (matrix, L x N)
%   nitot = plasma total ion density, cm^-3 
%   Ti = plasma ion temperature, keV 
%   Te = plasma electron temperature, keV 
%  Field ion(s):
%   Zf = plasma charge (matrix, L x N)
%   Af = plasma atomic number (matrix, L x N)
%  Varargin (controls): note, if one is used, all must be defined.
%   units = 'erg/cm', 'MeV/um', ['MeV/(mg/cm^2)'], 
%   TeDegeneracy = 0, [1];   if 1, correct electron temperature for electron degeneracy
%   diffStopping = [0], 1;   if 1, output stopping powers for all species separately
%   mimicFredrick = [0], 1;  if 1, mimic math in Fredrick's 
%                            stoppingpower9.exe program.  
%                            This has been benchmarked and gives basically
%                            identical results; note that there's a
%                            significant difference from not selecting this
%                            option.  Differences in: vf, vt, ur, xtf, cc, logLambda.
%Outputs:
%   stppwr = stopping power (units set by input, default = 'MeV/(mg/cm^2)')
%            (matrix size L x 1 (default, summed) or L x N+1 (output differential stopping)
%   By convention, stopping power output is positive.
%
%NOTES from Alex on appropriate thermal velocity (currently implemented):
%   1. In the Coulomb log through u, one should use 8kT/pi*m
%   2. In the Chandrasekhar function G through x^t/f, one should use 2kT/m
%   3. In the collective effects term, also through x^t/f, one should use kT/m

if isempty(varargin), 
    units = 'MeV/(mg/cm^2)';
    TeDegeneracy = 1; %correct electron temperature for electron degeneracy
    diffStopping = 0; %output stopping powers for different species
    mimicFredrick = 0; %use Fredrick's math
    quiet = 0;
elseif length(varargin)==4,
    units = varargin{1};
    TeDegeneracy = varargin{2};
    diffStopping = varargin{3};
    mimicFredrick = varargin{4};
    quiet = 1;
else 
    fprintf('Warning: incorrect controls declaration, using default.\n');
    units = 'MeV/(mg/cm^2)';
    TeDegeneracy = 1; %correct electron temperature for electron degeneracy
    diffStopping = 0; %output stopping powers for different species
    mimicFredrick = 0; %use Fredrick's math
    quiet = 0;
end

%% test values
%{
Et = [1,2,3,4,5];
L = length(Et);
Zt = 1*ones(L,1);
At = 1*ones(L,1);
nitot = 1e24*ones(L,1);
Ti = 1*ones(L,1);
Te = Ti;
Zf = repmat([6,1],[L,1]);
Af = repmat([12,1],[L,1]);
fi = repmat([1,1.4],[L,1]);

units = 'MeV/(mg/cm^2)';
diffStopping = 1;
TeDegeneracy = 0;
mimicFredrick = 0;
%}
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

%% make sure shape is correct:
L = length(Et);     %number of tests
N = size(Zf,2);     %number of plasma ion species
%Dimension 1: test scenarios
%Dimension 2: plasma ions
Et = reshape(Et,L,1); 

Zf = reshape(Zf,L,N); 
Zmat = cat(2,Zf,ones(L,1));
Af = reshape(Af,L,N); 
Afmat = cat(2,Af,me/mp*ones(L,1));

Atmat = repmat(At,[1,N+1]);
Ztmat = repmat(Zt,[1,N+1]);

%% Density, temperature
fi = fi./repmat(sum(fi,2),[1,N]); %renormalize fi
ni = repmat(nitot,[1,N]).*fi;
ne = sum(ni.*Zmat(:,1:N),2); %electron density, cm^-3
nmat = cat(2,ni,ne);

Ti = reshape(Ti,L,1);
Te = reshape(Te,L,1);

if TeDegeneracy,
    if mimicFredrick,
        Te = LiPetrassoTeffFredrick(ne,Te);
    else
        Te = LiPetrassoTeff(ne,Te,quiet);
    end
end

Tmat = cat(2,repmat(Ti,[1,N]),Te);

%% Calculate velocities, reduced quantities, Chandrasekhar function (xtf)
% values that are calculated separately for each species are matrices L x N
% everything else is vectors length L
% NOTE:  vf should be maxwellian average (8T/pi/m)^(1/2); xtf should use
% ratio of particle energy (i.e. "v_xtf" = 2T/m)
Ar = Atmat.*Afmat./(Atmat+Afmat); %reduced atm number

if ~mimicFredrick,
    vt = c*sqrt(((At.*mpc2+Et).^2-(At*mpc2).^2)./(At*mpc2+Et).^2); % cm/sec, velocity of test particle
    vf = sqrt(8/pi*kBkeV*Tmat./(Afmat*mp)); %sqrt(8/pi*kBkeV*Tmat./(Afmat*mp));  % cm/sec, average speed of field ion, electron
else
    vt = sqrt(2*Et./(At*mpc2))*c;
    vf = sqrt(3*(Tmat/1000)./(Afmat*mpc2))*c;
end
vtmat = repmat(vt,[1,N+1]);

if ~mimicFredrick,
    %Better approximation of relative velocity, <u> = <vt-vf>; makes assumption that vf^2 = <v>^2 = 8T/pi/m
    ur = vf/2.*exp(-(4/pi)*vtmat.^2./vf.^2)+vtmat.*(1+(pi/8)*vf.^2./vtmat.^2).*erf(sqrt(4/pi*vtmat.^2./vf.^2)); % cm/sec, relative velocity (average)
    xtf = vtmat.^2./(vf.^2*(2*pi/8)); % test particle/field particle energy ratio, uses vf = sqrt(2*T/m)
else
    %Basic approximation of relative velocity, <u> = <(vt-vf)^2>^(1/2); differs by at most 2.5%
    ur = sqrt(vtmat.^2 + vf.^2); % cm/sec, relative velocity (average)
    xtf = (Afmat./Atmat).*(repmat(Et,[1,N+1])./(Tmat/1000)); % test particle/field particle energy ratio
end

%% Calculate plasma parameters
if ~mimicFredrick, 
    lambdaDebye = (4*pi*e^2*sum(nmat.*Zmat.^2./(kBkeV*Tmat),2)).^(-1/2); %cm,  CLASSICAL
    lambdaDebyemat = repmat(lambdaDebye,[1,N+1]);

    pperp = Ztmat.*Zmat*e^2./(Ar*mp.*ur.^2); %cm
    pquantum = hbar./(2*Ar*mp.*ur);
    pmin = sqrt(pperp.^2 + pquantum.^2); %cm

    %logLambda = log(lambdaDebyemat./pmin);           %approximation > 1
    logLambda = 0.5*log(1+(lambdaDebyemat./pmin).^2); %Trubnikov, Eqn 1.5
else
    lambdaDebyemat = (4*pi*e^2*nmat.*Zmat.^2./(kBkeV*Tmat)).^(-1/2); %cm,  CLASSICAL
    logLambda = 0.5*log((1.06493e14*(Tmat/1000).*(Ar*mpc2).^2.*(ur/c).^2./(nmat/1e24))./(4*Ztmat.^2.*Zmat.^2./(ur/c).^2 + 137^2));
end

omegapf = sqrt(4*pi*nmat.*Zmat.^2*e^2./(Afmat*mp)); % sec, plasma frequency

%% Calculate G(xtf)
%mu is the Maxwell integral, int(2*e^(-x)*x^(1/2)/pi^(1/2),x,0,xtf), which
%has an evaluated form:
mu = erf(xtf.^(1/2)) - (2*xtf.^(1/2).*exp(-xtf))/pi^(1/2); % Maxwell integral from 0 to x
DmuDxtf = 2*exp(-xtf).*xtf.^(1/2)/pi^(1/2);

G = mu - Afmat./Atmat.*(DmuDxtf - (mu + DmuDxtf)./logLambda);

%% Calculate stopping power
% NOTE: use definition v_t^2 = T/m in collective-effects term; v_t^2 = 2T/m
% was used previously
xtf_collective = xtf*2; 

if ~mimicFredrick, 
    %After C.K. Li thesis, Appendix C eqn 8 (pg 209)
    logCollective = besselk(1,xtf_collective.^(-1/2)).*besselk(0,xtf_collective.^(-1/2))./xtf_collective.^(1/2);
else
    %Limit of above for xtf>>1
    cc = 1.23; %probably a typo???
    logCollective = (xtf_collective>=(cc^-2)).*log(cc*xtf_collective.^(1/2));
end
stppwr = Ztmat.^2*e^2./vtmat.^2.*omegapf.^2 ...
   .*(G.*logLambda + logCollective); % erg/cm

stppwr(stppwr<0) = 0;       %throw out bogus values at low energy
stppwr(fi==0) = 0;          %throw out values if ion fraction disappears
%% modify units and output
switch units
    case 'erg/cm'
        output = stppwr*1;
    case 'MeV/um'
        output = stppwr*(62.4150974);
    case 'MeV/(mg/cm^2)'
        output = stppwr*(6.24150974e5)... %MeV/cm
            ./repmat(sum(Afmat.*nmat*mp,2),[1,N+1])... %MeV/(g/cm^2)
            /1000;              %MeV/(mg/cm^2)
    otherwise
        error('Unknown units request.  Please input: ''erg/cm'',''MeV/um'', ''MeV/(mg/cm^2)''.\n');
end

if ~diffStopping,
    output = sum(output,2); %output only total stopping power
end

output(isnan(output)) = 0;
