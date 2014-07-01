function [sBER] = ber(snr)
a = round(rand(1,10000));

b=40;
register1=[1 1 1 1 1];
code1=zeros(1,b);
for i=1:b
   temp = mod(register1(2)+register1(5),2);
   code1(i) = 2*register1(5)-1;
   for j=5:-1:2
      register1(j)=register1(j-1);
   end
   register1(1) = temp;
end
m_sequence_1=code1;
m_sequence_1=m_sequence_1';
m_sequence_1 =m_sequence_1>0;
l=length(m_sequence_1);
czas=0:1:(l-1);
stem(czas,m_sequence_1)
m_sequence_1=num2str(m_sequence_1);
l = length(a);
b = a -'0';
c=m_sequence_1-'0';
c=c';
k=1;
g=length(m_sequence_1);
for i=1:l
for j=1:g
    
    spread(1,k)=xor(b(1,i),c(1,j));
    k=k+1;
end
end

l=length(spread);
z=0:1:(l-1);


NRZ = 2*a-1;                                        % Kodowanie NRZ
T = 1;                                              % Czas trwania bitu (okres)
fc = 3/T;                                           % Czestotlowosc nosna
                                         % Czestotlowosc nosna
t = linspace(0,length(NRZ),length(NRZ)*100);      % probki, czasy
N = length(t);                                      % Liczba probek
Lpnb = N/length(NRZ);                              % Liczba probek na bit
                    % powtarzamy bity Lpnb razy
dod_dane2 = repmat(NRZ',1,Lpnb);
                                   % Transpozycja wierszy i kolumn
d_d2 = dod_dane2';
                                    % Konwersja do jednego wiersza (100x1, 100x0, itd..)
d_d2 = d_d2(:)';     
nosna = sin(pi*fc*t);                               % Fala nosna 
bpsk_mod = d_d2.*nosna;    


szum = awgn(bpsk_mod,snr);    

odfiltr=szum.*nosna;
for i=1:length(odfiltr)
    if odfiltr(i)>0
        bpsk_demod(i)=1;
    else
        bpsk_demod(i)=0;
    end
end
% Bledy - zwraca ich ilosc
y = bpsk_demod;                  
w = real(y)>0;
z=zeros(1,length(NRZ));
j=z;
bit=z;
for q=0:length(NRZ)-1
    for b=1:15
        if w(100*q+6*b)==0
            z(q+1)=z(q+1)+1;
        else j(q+1)=j(q+1)+1;
        end     
    end
    
end
for q=0:length(NRZ)-1
    for b=1:15
        if z(q+1)>j(q+1)
            wnew(15*q+b)=0;
            bit(q+1)=0;
        else
            wnew(15*q+b)=1;
            bit(q+1)=1;
        end
    end
end

nErr = size(find([a- bit]),2);
sBER = nErr/length(a)

function y=awgn(varargin)
%AWGN Add white Gaussian noise to a signal.
%   Y = AWGN(X,SNR) adds white Gaussian noise to X.  The SNR is in dB.
%   The power of X is assumed to be 0 dBW.  If X is complex, then 
%   AWGN adds complex noise.
%
%   Y = AWGN(X,SNR,SIGPOWER) when SIGPOWER is numeric, it represents 
%   the signal power in dBW. When SIGPOWER is 'measured', AWGN measures
%   the signal power before adding noise.
%
%   Y = AWGN(X,SNR,SIGPOWER,STATE) resets the state of RANDN to STATE.
%
%   Y = AWGN(..., POWERTYPE) specifies the units of SNR and SIGPOWER.
%   POWERTYPE can be 'db' or 'linear'.  If POWERTYPE is 'db', then SNR
%   is measured in dB and SIGPOWER is measured in dBW.  If POWERTYPE is
%   'linear', then SNR is measured as a ratio and SIGPOWER is measured
%   in Watts.
%
%   Example: To specify the power of X to be 0 dBW and add noise to produce
%            an SNR of 10dB, use:
%            X = sqrt(2)*sin(0:pi/8:6*pi);
%            Y = AWGN(X,10,0);
%
%   Example: To specify the power of X to be 0 dBW, set RANDN to the 1234th
%            state and add noise to produce an SNR of 10dB, use:
%            X = sqrt(2)*sin(0:pi/8:6*pi);
%            Y = AWGN(X,10,0,1234);
%
%   Example: To specify the power of X to be 3 Watts and add noise to
%            produce a linear SNR of 4, use:
%            X = sqrt(2)*sin(0:pi/8:6*pi);
%            Y = AWGN(X,4,3,'linear');
%
%   Example: To cause AWGN to measure the power of X, set RANDN to the 
%            1234th state and add noise to produce a linear SNR of 4, use:
%            X = sqrt(2)*sin(0:pi/8:6*pi);
%            Y = AWGN(X,4,'measured',1234,'linear');
%
%   See also WGN, RANDN, and BSC.

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.9.4.2 $  $Date: 2004/04/12 23:00:22 $ 

% --- Initial checks
error(nargchk(2,5,nargin));

% --- Value set indicators (used for the string flags)
pModeSet    = 0;
measModeSet = 0;

% --- Set default values
reqSNR   = [];
sig      = [];
sigPower = 0;
pMode    = 'db';
measMode = 'specify';
state    = [];

% --- Placeholder for the signature string
sigStr = '';

% --- Identify string and numeric arguments
for n=1:nargin
   if(n>1)
      sigStr(size(sigStr,2)+1) = '/';
   end
   % --- Assign the string and numeric flags
   if(ischar(varargin{n}))
      sigStr(size(sigStr,2)+1) = 's';
   elseif(isnumeric(varargin{n}))
      sigStr(size(sigStr,2)+1) = 'n';
   else
      error('Only string and numeric arguments are allowed.');
   end
end

% --- Identify parameter signatures and assign values to variables
switch sigStr
   % --- awgn(x, snr)
   case 'n/n'
      sig      = varargin{1};
      reqSNR   = varargin{2};

   % --- awgn(x, snr, sigPower)
   case 'n/n/n'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      sigPower = varargin{3};

   % --- awgn(x, snr, 'measured')
   case 'n/n/s'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      measMode = lower(varargin{3});

      measModeSet = 1;

   % --- awgn(x, snr, sigPower, state)
   case 'n/n/n/n'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      sigPower = varargin{3};
      state    = varargin{4};

   % --- awgn(x, snr, 'measured', state)
   case 'n/n/s/n'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      measMode = lower(varargin{3});
      state    = varargin{4};

      measModeSet = 1;

   % --- awgn(x, snr, sigPower, 'db|linear')
   case 'n/n/n/s'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      sigPower = varargin{3};
      pMode    = lower(varargin{4});

      pModeSet = 1;

   % --- awgn(x, snr, 'measured', 'db|linear')
   case 'n/n/s/s'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      measMode = lower(varargin{3});
      pMode    = lower(varargin{4});

      measModeSet = 1;
      pModeSet    = 1;

   % --- awgn(x, snr, sigPower, state, 'db|linear')
   case 'n/n/n/n/s'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      sigPower = varargin{3};
      state    = varargin{4};
      pMode    = lower(varargin{5});

      pModeSet = 1;

   % --- awgn(x, snr, 'measured', state, 'db|linear')
   case 'n/n/s/n/s'
      sig      = varargin{1};
      reqSNR   = varargin{2};
      measMode = lower(varargin{3});
      state    = varargin{4};
      pMode    = lower(varargin{5});

      measModeSet = 1;
      pModeSet    = 1;

   otherwise
      error('Syntax error.');
end   

% --- Parameters have all been set, either to their defaults or by the values passed in,
%     so perform range and type checks

% --- sig
if(isempty(sig))
   error('An input signal must be given.');
end

if(ndims(sig)>2)
   error('The input signal must have 2 or fewer dimensions.');
end

% --- measMode
if(measModeSet)
   if(~strcmp(measMode,'measured'))
      error('The signal power parameter must be numeric or ''measured''.');
   end
end

% --- pMode
if(pModeSet)
   switch pMode
   case {'db' 'linear'}
   otherwise
      error('The signal power mode must be ''db'' or ''linear''.');
   end
end

% -- reqSNR
if(any([~isreal(reqSNR) (length(reqSNR)>1) (length(reqSNR)==0)]))
   error('The signal-to-noise ratio must be a real scalar.');
end

if(strcmp(pMode,'linear'))
   if(reqSNR<=0)
      error('In linear mode, the signal-to-noise ratio must be > 0.');
   end
end

% --- sigPower
if(~strcmp(measMode,'measured'))

   % --- If measMode is not 'measured', then the signal power must be specified
   if(any([~isreal(sigPower) (length(sigPower)>1) (length(sigPower)==0)]))
      error('The signal power value must be a real scalar.');
   end
   
   if(strcmp(pMode,'linear'))
      if(sigPower<0)
         error('In linear mode, the signal power must be >= 0.');
      end
   end

end

% --- state
if(~isempty(state))
   if(any([~isreal(state) (length(state)>1) (length(state)==0) any((state-floor(state))~=0)]))
      error('The State must be a real, integer scalar.');
   end
end

% --- All parameters are valid, so no extra checking is required

% --- Check the signal power.  This needs to consider power measurements on matrices
if(strcmp(measMode,'measured'))
   sigPower = sum(abs(sig(:)).^2)/length(sig(:));

   if(strcmp(pMode,'db'))
      sigPower = 10*log10(sigPower);
   end
end

% --- Compute the required noise power
switch lower(pMode)
   case 'linear'
      noisePower = sigPower/reqSNR;
   case 'db'
      noisePower = sigPower-reqSNR;
      pMode = 'dbw';
end

% --- Add the noise
if(isreal(sig))
   opType = 'real';
else
   opType = 'complex';
end

y = sig+wgn(size(sig,1), size(sig,2), noisePower, 1, state, pMode, opType);

function y = wgn(varargin)
%WGN Generate white Gaussian noise.
%   Y = WGN(M,N,P) generates an M-by-N matrix of white Gaussian noise.
%   P specifies the power of the output noise in dBW.
%
%   Y = WGN(M,N,P,IMP) specifies the load impedance in Ohms.
%
%   Y = WGN(M,N,P,IMP,STATE) resets the state of RANDN to STATE.
%
%   Additional flags that can follow the numeric arguments are:
%
%   Y = WGN(..., POWERTYPE) specifies the units of P.  POWERTYPE can
%   be 'dBW', 'dBm' or 'linear'.  Linear power is in Watts.
%
%   Y = WGN(..., OUTPUTTYPE); Specifies the output type.  OUTPUTTYPE can
%   be 'real' or 'complex'.  If the output type is complex, then P
%   is divided equally between the real and imaginary components.
%
%   Example: To generate a 1024-by-1 vector of complex noise with power
%            of 5 dBm across a 50 Ohm load, use:
%            Y = WGN(1024, 1, 5, 50, 'dBm', 'complex');
%
%   Example: To generate a 256-by-5 matrix of real noise with power
%            of 10 dBW across a 1 Ohm load, use:
%            Y = WGN(256, 5, 10, 'real');
%
%   Example: To generate a 1-by-10 vector of complex noise with power
%            of 3 Watts across a 75 Ohm load, use:
%            Y = WGN(1, 10, 3, 75, 'linear', 'complex');
%
%   See also RANDN, AWGN.

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.11.4.3 $  $Date: 2004/04/12 23:01:28 $

% --- Initial checks
error(nargchk(3,7,nargin));

% --- Value set indicators (used for the strings)
pModeSet    = 0;
cplxModeSet = 0;

% --- Set default values
p        = [];
row      = [];
col      = [];
pMode    = 'dbw';
imp      = 1;
cplxMode = 'real';
seed     = [];

% --- Placeholders for the numeric and string index values
numArg = [];
strArg = [];

% --- Identify string and numeric arguments
%     An empty in position 4 (Impedance) or 5 (Seed) are considered numeric
for n=1:nargin
   if(isempty(varargin{n}))
      switch n
      case 4
         if(ischar(varargin{n}))
            error('The default impedance should be marked by [].');
         end;
         varargin{n} = imp; % Impedance has a default value
      case 5
         if(ischar(varargin{n}))
            error('The default seed should be marked by [].');
         end;
         varargin{n} = [];  % Seed has no default
      otherwise
         varargin{n} = '';
      end;
   end;

   % --- Assign the string and numeric vectors
   if(ischar(varargin{n}))
      strArg(size(strArg,2)+1) = n;
   elseif(isnumeric(varargin{n}))
      numArg(size(numArg,2)+1) = n;
   else
      error('Only string and numeric arguments are allowed.');
   end;
end;

% --- Build the numeric argument set
switch(length(numArg))

   case 3
      % --- row is first (element 1), col (element 2), p (element 3)

      if(all(numArg == [1 2 3]))
         row    = varargin{numArg(1)};
         col    = varargin{numArg(2)};
         p      = varargin{numArg(3)};
      else
         error('Illegal syntax.')
      end;

   case 4
      % --- row is first (element 1), col (element 2), p (element 3), imp (element 4)
      %

      if(all(numArg(1:3) == [1 2 3]))
         row    = varargin{numArg(1)};
         col    = varargin{numArg(2)};
         p      = varargin{numArg(3)};
         imp    = varargin{numArg(4)};
      else
         error('Illegal syntax.')
      end;

   case 5
      % --- row is first (element 1), col (element 2), p (element 3), imp (element 4), seed (element 5)

      if(all(numArg(1:3) == [1 2 3]))
         row    = varargin{numArg(1)};
         col    = varargin{numArg(2)};
         p      = varargin{numArg(3)};
         imp    = varargin{numArg(4)};
         seed   = varargin{numArg(5)};
      else
         error('Illegal syntax.');
      end;
   otherwise
      error('Illegal syntax.');
end;

% --- Build the string argument set
for n=1:length(strArg)
   switch lower(varargin{strArg(n)})
   case {'dbw' 'dbm' 'linear'}
      if(~pModeSet)
         pModeSet = 1;
         pMode = lower(varargin{strArg(n)});
      else
         error('The Power mode must only be set once.');
      end;
   case {'db'}
      error('Incorrect power mode passed in.  Please use ''dBW'', ''dBm'', or ''linear.''');
   case {'real' 'complex'}
      if(~cplxModeSet)
         cplxModeSet = 1;
         cplxMode = lower(varargin{strArg(n)});
      else
         error('The complexity mode must only be set once.');
      end;
   otherwise
      error('Unknown option passed in.');
   end;
end;

% --- Arguments and defaults have all been set, either to their defaults or by the values passed in
%     so, perform range and type checks

% --- p
if(isempty(p))
   error('The power value must be a real scalar.');
end;

if(any([~isreal(p) (length(p)>1) (length(p)==0)]))
   error('The power value must be a real scalar.');
end;

if(strcmp(pMode,'linear'))
   if(p<0)
      error('In linear mode, the required noise power must be >= 0.');
   end;
end;

% --- Dimensions
if(any([isempty(row) isempty(col) ~isscalar(row) ~isscalar(col)]))
   error('The required dimensions must be real, integer scalars > 1.');
end;

if(any([(row<=0) (col<=0) ~isreal(row) ~isreal(col) ((row-floor(row))~=0) ((col-floor(col))~=0)]))
   error('The required dimensions must be real, integer scalars > 1.');
end;

% --- Impedance
if(any([~isreal(imp) (length(imp)>1) (length(imp)==0) any(imp<=0)]))
   error('The Impedance value must be a real scalar > 0.');
end;

% --- Seed
if(~isempty(seed))
   if(any([~isreal(seed) (length(seed)>1) (length(seed)==0) any((seed-floor(seed))~=0)]))
      error('The State must be a real, integer scalar.');
   end;
end;

% --- All parameters are valid, so no extra checking is required
switch lower(pMode)
   case 'linear'
      noisePower = p;
   case 'dbw'
      noisePower = 10^(p/10);
   case 'dbm'
      noisePower = 10^((p-30)/10);
end;

% --- Generate the noise
if(~isempty(seed))
   randn('state',seed);
end;

if(strcmp(cplxMode,'complex'))
   z = randn(2*row,col);
   y = (sqrt(imp*noisePower/2))*(z(1:row,:)+j*z(row+1:end,:));
else
   y = (sqrt(imp*noisePower))*randn(row,col);
end;

