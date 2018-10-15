function [causality,coherence,power] = compute_parCausality(X,fs,freq,p)
%Usage:  [causality,coherence,power] = compute_parCausality(X,fs,freq,p); 
%-------------------------------------------------------------------------------------------------
%Inputs: X = bivariate or multivariate signals in the form of 3D matrix  for time. trial. channel.  
%               fs = data sampling rate in Hz
%               freq = a vector of frequencies to compute causality spectra (e.g. 0:1:100) 
%               p = model order for autoregressive fitting
%Outputs: causality = Granger causality between all pairs of signals (channels), e.g. 
%               :             causality(:,2,1) means causality from 2 to 1
%               :             auto-causality are set to zero, i.e. causality(:,k,k) = 0 for all k
%               : coherence = coherence between all pairs of signals (form: frequency. channel. channel)
%               : power      = 1-sided auto-spectra (form: frequency. channel)
%-------------------------------------------------------------------------------------------------
%This function utilizes armorf.m, spectrum.m written by Y. Chen & colleagues
%   Ref: [1]. M. Morf, etal, Recursive Multichannel Maximum Entropy Spectral Estimation,
%                    IEEE trans. GeoSci. Elec., 1978, Vol.GE-16, No.2, pp85-94.
%                    S. Haykin, Nonlinear Methods of Spectral Analysis, 2nd Ed.
%                    Springer-Verlag, 1983, Chapter 2
%            [2]. Mingzhou Ding, et al. Biol. Cybern. 83, 2000:35-45.
%            [3]. Brovelli, et. al., PNAS 101, 9849-9854 (2004).
%-------------------------------------------------------------------------------------------------
%Flow of computation:  from signals (X) to 
%                                                1. parameters (A, Z) by model fitting 
%                                                2. transfer function (H), auto-& cross spectra (S)
%                                                3. Granger causality, coherence, 1-sided power (2*auto-spectra) 
% M. Dhamala, Univ. of Florida, Gainesville, Aug 2006.
%-------------------------------------------------------------------------------------------------

[Nt, Ntr,Nc] = size(X); %Nt = number of timepoints, Ntr = trials, Nc = channels 
 
x = reshape(X,Nt*Ntr,Nc); x = x'; %' x in the form of channel. (time x trials)

[A, Z]=armorf(x,Ntr,Nt,p); %parameters by autoregressive fitting

[S,H] = AZ2spectra(A,Z,p,freq,fs); 

if nargout>1
   spectra = permute(S,[3 1 2]);
   coherence = S2coh(spectra);
  if nargout>2
       for ichan = 1: Nc,
             power(:,ichan) = 2*spectra(:,ichan,ichan);%1-sided power
       end
  end
end

if size(X,3)<3, % for pairwise causality of 2-channels, use H, S, Z already computed
     causality = hz2causality(H,S,Z,fs); 
else % for more than 2 channels, fit 2-channels at a time and find pairwise causality 
     for ii = 1: Nc-1,
         for jj = ii+1: Nc, 
               [A2,Z2] = armorf(x([ii jj],:),Ntr,Nt,p); [S2,H2] = AZ2spectra(A2,Z2,p,freq,fs);
               cs = hz2causality(H2,S2,Z2,fs); 
               causality(:,ii,jj) = cs(:,1,2); causality(:,jj,ii) = cs(:,2,1);
         end
         causality(:,ii,ii) = 0; % self-causality is set to zero
     end 
     causality(:,Nc,Nc) = 0; %self-causality of the last channel is set to zero too   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------------------
% Functions below : armorf.m, AZ2spectra.m, spectrum.m, S2coh, hz2causality.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-------------------------------------------------------

%-------------------------------armorf.m-------------------------------------------------------

function varargout = armorf(x,Nr,Nl,p);
%ARMORF   AR parameter estimation via LWR method by Morf modified.
%   x is a matrix whose every row is one variable's time series
%   Nr is the number of realizations, Nl is the length of every realization
%   If the time series are stationary long, just let Nr=1, Nl=length(x)
%   p is the order of AR model
%
%   A = ARMORF(X,NR,NL,P) returns the polynomial coefficients A corresponding to 
%     the AR model estimate of matrix X using Morf's method.
%
%   [A,E] = ARMORF(...) returns the final prediction error E (the
%   covariance matrix of the white noise of the AR model).
%
%   [A,E,K] = ARMORF(...) returns the vector K of reflection 
%     coefficients (parcor coefficients).
%
%   Ref: M. Morf, etal, Recursive Multichannel Maximum Entropy Spectral Estimation,
%              IEEE trans. GeoSci. Elec., 1978, Vol.GE-16, No.2, pp85-94.
%        S. Haykin, Nonlinear Methods of Spectral Analysis, 2nd Ed.
%              Springer-Verlag, 1983, Chapter 2

% Initialization
[L,N]=size(x);
R0=zeros(L,L);
R0f=R0;
R0b=R0;
pf=R0;
pb=R0;
pfb=R0;
ap(:,:,1)=R0;
bp(:,:,1)=R0;
En=R0;
for i=1:Nr
    En=En+x(:,(i-1)*Nl+1:i*Nl)*x(:,(i-1)*Nl+1:i*Nl)';
    ap(:,:,1)=ap(:,:,1)+x(:,(i-1)*Nl+2:i*Nl)*x(:,(i-1)*Nl+2:i*Nl)';        
    bp(:,:,1)=bp(:,:,1)+x(:,(i-1)*Nl+1:i*Nl-1)*x(:,(i-1)*Nl+1:i*Nl-1)';
end
ap(:,:,1) = inv((chol(ap(:,:,1)/Nr*(Nl-1)))');
bp(:,:,1) = inv((chol(bp(:,:,1)/Nr*(Nl-1)))');
for i=1:Nr
    efp = ap(:,:,1)*x(:,(i-1)*Nl+2:i*Nl);
    ebp = bp(:,:,1)*x(:,(i-1)*Nl+1:i*Nl-1);
    pf = pf + efp*efp';
    pb = pb + ebp*ebp';
    pfb = pfb + efp*ebp';
end
En = chol(En/N)'; % Covariance of the noise

% Initial output variables
coeff = [];%  Coefficient matrices of the AR model
kr=[];  % reflection coefficients

for m=1:p
   % Calculate the next order reflection (parcor) coefficient
   ck = inv((chol(pf))')*pfb*inv(chol(pb));
   kr=[kr,ck];
   % Update the forward and backward prediction errors
   ef = eye(L)- ck*ck';
   eb = eye(L)- ck'*ck;
     
   % Update the prediction error
   En = En*chol(ef)';
   E = (ef+eb)./2;   
   
   % Update the coefficients of the forward and backward prediction errors
   ap(:,:,m+1) = zeros(L);
   bp(:,:,m+1) = zeros(L);
   pf = zeros(L);
   pb = zeros(L);
   pfb = zeros(L);

   for i=1:m+1       
       a(:,:,i) = inv((chol(ef))')*(ap(:,:,i)-ck*bp(:,:,m+2-i));
       b(:,:,i) = inv((chol(eb))')*(bp(:,:,i)-ck'*ap(:,:,m+2-i));
   end
   for k=1:Nr
       efp = zeros(L,Nl-m-1);
       ebp = zeros(L,Nl-m-1);
       for i=1:m+1
           k1=m+2-i+(k-1)*Nl+1;
           k2=Nl-i+1+(k-1)*Nl;
           efp = efp+a(:,:,i)*x(:,k1:k2);
           ebp = ebp+b(:,:,m+2-i)*x(:,k1-1:k2-1);
       end
       pf = pf + efp*efp';
       pb = pb + ebp*ebp';
       pfb = pfb + efp*ebp';
   end
   ap = a;
   bp = b;
end
for j=1:p
    coeff = [coeff,inv(a(:,:,1))*a(:,:,j+1)];
end

varargout{1} = coeff;
if nargout >= 2
    varargout{2} = En*En';
end
if nargout >= 3
    varargout{3} = kr;
end

%----------------------------- AZ2spectra.m --------------------------------    

function [S,H] = AZ2spectra(A,Z,p,freq,fs);
f_ind = 0;
for f = freq,
      f_ind = f_ind+1;
      [Stmp,Htmp] = spectrum(A,Z,p,f,fs);       
      H(:,:,f_ind) = Htmp;
      S(:,:,f_ind) = Stmp; %auto-& cross-spectra 
end
for ind = 1:size(Z,1),
      S(ind,ind,:) = real(S(ind,ind,:)); %avoiding numerical errors
end

%------------------------------spectrum.m-----------------------------------

function [S,H] = spectrum(A,Z,M,f,fs);
% Get the coherence spectrum
N = size(Z,1);
H = eye(N,N); % identity matrix
for m = 1 : M
    H = H + A(:,(m-1)*N+1:m*N)*exp(-i*m*2*pi*f/fs);
    % Multiply f in the exponent by sampling interval (=1/fs). See Shiavi
end
H = inv(H);
S = H*Z*H'/fs;
% One has to multiply HZH' by sampling interval (=1/fs)
%to get the properly normalized spectral density. See Shiavi.
%To get 1-sided power spectrum, multiply S by 2.

%--------------------------------------------S2coh.m -----------------------

function coh = S2coh(S); 
%Input: S auto-& cross pectra in the form: frequency. channel. channel 
%Output: coh (Coherence) in the form: frequency. channel. channel 
%M. Dhamala, UF, August 2006.

Nc = size(S,2);
for ii = 1: Nc,
   for jj = 1: Nc,
       coh(:,ii,jj) = real(abs(S(:,ii,jj)).^2./(S(:,ii,ii).*S(:,jj,jj)));
   end
end

%-------------------------------hz2causality.m ---------------------------- 

function causality = hz2causality(H,S,Z,fs);
%Usage: causality = hz2causality(H,S,Z,fs);
%Inputs: H = transfer function, S = 3-D spectral matrix;
%        Z = noise covariance,  fs = sampling rate
%Outputs: causality (Granger causality between all channels)
%               : auto-causality spectra are set to zero
% Reference: Brovelli, et. al., PNAS 101, 9849-9854 (2004).
%M. Dhamala, UF, August 2006.

Nc = size(H,2);

for ii = 1: Nc,
    for jj = 1: Nc,
          if ii ~=jj,
              zc = Z(jj,jj) - Z(ii,jj)^2/Z(ii,ii);
              numer = abs(S(ii,ii,:));
              denom = abs(S(ii,ii,:)-zc*abs(H(ii,jj,:)).^2/fs);
              causality(jj,ii,:) = log(numer./denom);
          end
    end
    causality(ii,ii,:) = 0;%self-causality set to zero
end
causality = permute(causality,[3 1 2]); %freq x channel from x channel to
