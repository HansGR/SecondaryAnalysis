clear all;
i = 1; 
[Y1(i), Y1unc(i), Y2n(i), Y2nunc(i), Y2p(i), Y2punc(i), Ti(i), Tiunc(i), rho(i), rhounc(i)] = deal(...
     3.01e11, 1.13e10, 8.57e7, 1.0e7, 2.01e8, 0.4e8, 5.4, 0.2, 0.356, 0.132);   %, N110131

i = 2; 
[Y1(i), Y1unc(i), Y2n(i), Y2nunc(i), Y2p(i), Y2punc(i), Ti(i), Tiunc(i), rho(i), rhounc(i)] = deal(...
     2.52e11, 0.184e11, 7.6e7, 2.8e7, 1.7e8, 0.3e8, 3.9, 0.3, 0.403, 0.126);    %, N130129

%%% NIF exploding pushers:
%     3.01e11, 1.13e10, 8.57e7, 1.0e7, 2.01e8, 0.4e8, 5.4, 0.2, 0.356, 0.132);   %, N110131
%     2.52e11, 0.184e11, 7.6e7, 2.8e7, 1.7e8, 0.3e8, 3.9, 0.3, 0.403, 0.126);    %, N130129

Te = [0.67,0.75].*Ti;
Zf = repmat([1,14,8],[length(Ti),1]);
Af = repmat([2,28,16],[length(Ti),1]);
fi = repmat([1,0,0],[length(Ti),1]);
HSorUni = 1; 

%a = FindRhoRFromSecondaryYield_old(Y1,Y1unc,Y2n, Y2nunc, Y2p, Y2punc, Ti, Te, rho, Zf, Af, fi, HSorUni, 1);
b = FindRhoRFromSecondaryYield(Y1,Y1unc,Y2n, Y2nunc, Y2p, Y2punc, Ti, Te, rho, Zf, Af, fi, HSorUni);


%% test Mix
clear all
[Y1, Y1unc, Y2n, Y2nunc, Y2p, Y2punc, Ti, Ti_unc, rho_0, rhounc_0] = deal(...
     3.01e11, 1.13e10, 8.57e7, 1.0e7, 2.01e8, 0.4e8, 5.4, 0.2, 0.356, 0.132);   %, N110131
%   2.52e11, 0.184e11, 7.6e7, 2.8e7, 1.7e8, 0.3e8, 3.9, 0.3, 0.403, 0.126);    %, N130129


Te = 0.67*Ti;
%Te = 0.75*Ti; 
Z0 = [1];
A0 = [2];
f0 = [1];
    f0 = f0/sum(f0);
Zmix = [14,8];
Amix = [28,16];
fmix = [1,2]; 
    fmix = fmix/sum(fmix);

mixfrac = [0.0:0.01:1]*0.2; %fraction of mix by atom, n_mix = n_fuel * mixfrac
mixfrac = mixfrac';
N = length(mixfrac);

fi = [repmat(f0,[N,1]),repmat(fmix,[N,1]).*repmat(mixfrac,[1,size(fmix,2)])];
HSorUni = 1; 
TeCorrection = 0;

a = FindMixFromSecondaries(Y1, Y1unc, Y2n, Y2nunc, Y2p, Y2punc, ...
    Ti, Te, rho_0, Z0, A0, Zmix, Amix, fi, HSorUni, 'TeCorrection', TeCorrection,'makeFigures',1);

%% test Te
clear all
[Y1, Y1unc, Y2n, Y2nunc, Y2p, Y2punc, Ti, Ti_unc, rho_0, rhounc_0] = deal(...
     22e12, 2e12, 13.5e10, 1.3e10, 14e9, 7e9, 3.6, 0.2, 2.96, 0.89);  %N130813 2-shock HDC


Te = 0.5:0.05:2.0;
Z0 = [1,2];
A0 = [2,3];
f0 = [1,0];

HSorUni = 1; 
TeCorrection = 0;

a = FindTeFromSecondaries(Y1, Y1unc, Y2n, Y2nunc, Y2p, Y2punc, ...
    Ti, rho_0, Z0, A0, f0, Te, HSorUni, 'TeCorrection', TeCorrection,'makeFigures',1);
