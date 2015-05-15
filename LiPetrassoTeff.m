function Teout = LiPetrassoTeff(ne, Te, varargin)
%Calculate the effective electron temperature due to degeneracy; 
%input:  (ne (cm^-3), Te (keV), [quiet = 0 [1], suppress messages])

%Commentary:
%{
Paul Drake talks about Fermi degenerate electrons in section 3.1.3 of his
book.  Eqn. 3.22 states that the electron pressure in a Fermi-degenerate
plasma is Pe = (2/3) n_e k_B T_e (F(3/2,alpha)/F(1/2,alpha)), where 
alpha = chemical potential mu/(k_B T_e).  

Damien's code for the degenerate electron temperature correction absorbs
all additional factors into T_e such that Pe = n_e k_B T_e_effective.
This code does the same, but with improved math.

Damien's approach was as follows:
1. find the classical value of F(1/2,alpha), based on ne, Te
2. find the value of alpha that produces this F(1/2,alpha)
3. use this alpha to calculate the Te correction.

This is conceptually problematic in that you are assuming the classical
limit at the beginning to calculate a non-classical phenomenon.  For highly
degenerate plasmas, this would underestimate the degeneracy effect.  

Drake includes a fit for alpha in terms of Theta = Te/Tdegenerate, from
Atzeni and Meyer-ter-Vehn due to Ichimaru (eqn 3.20).  This is used in this
code, which follows this procedure:
1. find Theta = Te/Td, from Te, ne
2. calculate Alpha using eqn 3.20
3. calculate Te-correction using F(1/2,alpha), F(3/2,alpha)
%}

if ~isempty(varargin),
    quiet = varargin{1};
else
    quiet = 1;
end

%%%%%%%%%%%%%% Constants
me = 9.10938291*10^-28; %g
mec2 = 0.510998903;     %MeV, electron rest mass
mp = 1.67262178*10^-24; %g
mpc2 = 938.272046;      %MeV, proton rest mass
h = 4.1356668e-12;      %MeV*ns, planck constant
c = 29.9792458;         %cm/ns,  speed of light
%%%%%%%%%%%%%%

%find Theta = Te/Tdegenerate, from Te, ne;  Drake eqn 3.16
Theta = (Te/1000).*(8*pi/3./ne).^(2/3).*(2*mec2/(h^2*c^2));

%find Alpha from Theta, using Ichimaru fit;  Drake eqn 3.20
alpha = -3/2*log(Theta) + log(4/3/sqrt(pi)) + ...
    (0.25054*Theta.^(-1.858)+0.072.*Theta.^(-1.858/2))./...
    (1+0.25054.*Theta.^(-0.858));

%calculate corrected electron temperature
%Te = Te_classical to within 0.1% for alpha <= -10
classical = (alpha<=-10);
Teout = zeros(size(Te));
Teout(classical) = Te(classical);

if (sum(~classical)>0)&&(~quiet), 
    fprintf('Note: non-classical plasma detected, electron degeneracy correction applied.\n'); 
end
Teout(~classical) = Te(~classical).*(2/3).*...
    F(3/2,alpha(~classical))./F(1/2,alpha(~classical));

%degenerate = (alpha>10);
%Teout(degenerate) = Te(degenerate)*(2/3).*F(3/2,1./(Theta(degenerate)))./F(1/2,1./(Theta(degenerate)));

function fout = F(n,alpha)
%If statement protects the integrand from approaching numerically small values
fout = zeros(size(alpha));
for i = 1:length(alpha),
    if alpha<0,
        integrand = @(x)(x.^n.*(exp(x)+exp(alpha(i))).^(-1)); %integrand, divided by exp(alpha(i))
        peak = max(alpha(i),2*n);
        fout(i) = quad(integrand,0,peak*10).*exp(alpha(i));
    else
        integrand = @(x)(x.^n.*(exp(x-alpha(i))+1).^(-1)); %integrand
        peak = max(alpha(i),2*n);
        fout(i) = quad(integrand,0,peak*10);
    end
end
