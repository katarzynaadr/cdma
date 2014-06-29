function varargout = cdmagui(varargin)
%CDMAGUI M-file for cdmagui.fig
%      CDMAGUI, by itself, creates a new CDMAGUI or raises the existing
%      singleton*.
%
%      H = CDMAGUI returns the handle to a new CDMAGUI or the handle to
%      the existing singleton*.
%
%      CDMAGUI('Property','Value',...) creates a new CDMAGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to cdmagui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CDMAGUI('CALLBACK') and CDMAGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CDMAGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cdmagui

% Last Modified by GUIDE v2.5 29-Jun-2014 13:35:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cdmagui_OpeningFcn, ...
                   'gui_OutputFcn',  @cdmagui_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cdmagui is made visible.
function cdmagui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for cdmagui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cdmagui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cdmagui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in p1.
function p1_Callback(hObject, eventdata, handles)
% hObject    handle to p1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zmienna_pomocnicza
get(handles.syg,'String')
textString=get(handles.syg,'String');
if(strcmp(textString,'0.')==1)&(zmienna_pomocnicza == 0)
    set(handles.syg,'String','1');
else
    textString=strcat(textString,'1');
    set(handles.syg,'String',textString)
end
zmienna_pomocnicza = 0;
guidata(hObject, handles);

% --- Executes on button press in p0.
function p0_Callback(hObject, eventdata, handles)
% hObject    handle to p0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zmienna_pomocnicza
get(handles.syg,'String')
textString=get(handles.syg,'String');
if(strcmp(textString,'0.')==1)&(zmienna_pomocnicza == 0)
    set(handles.syg,'String','0');
else
    textString=strcat(textString,'0');
    set(handles.syg,'String',textString)
end
zmienna_pomocnicza = 0;
guidata(hObject, handles);

% --- Executes on button press in przebieg.
function przebieg_Callback(hObject, eventdata, handles)
% hObject    handle to przebieg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a =get(handles.syg,'String');
guidata(hObject, handles);
l = length(a);
b =a-'0';
t=0:1:(l-1);
%b=rectpulse(a(1,:),l)
hold on
axes(handles.axes1)
stem(t,b)

% --- Executes on button press in sekwencja.
function sekwencja_Callback(hObject, eventdata, handles)
% hObject    handle to sekwencja (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%m-sequence (52)
register1=[1 1 1 1 1];
code1=zeros(1,31);
for i=1:31
   temp = mod(register1(2)+register1(5),2);
   code1(i) = 2*register1(5)-1;
   for j=5:-1:2
      register1(j)=register1(j-1);
   end
   register1(1) = temp;
end
%m-sequence (5432)
%register2=[1 1 1 1 1];%initial fill
%code2=zeros(1,31);
%for i=1:31
%   temp = mod(register2(2)+register2(3)+register2(4)+register2(5),2);
%   code2(i) = 2*register2(5)-1;
%   for j=5:-1:2
%      register2(j)=register2(j-1);
%   end
%   register2(1) = temp;
%end
m_sequence_1=code1;
%m_sequence_2=code2;
m_sequence_1=m_sequence_1';
%m_sequence_2=m_sequence_2';
m_sequence_1 =m_sequence_1>0;
%m_sequence_2 =m_sequence_2>0;
axes(handles.axes2)
stem(m_sequence_1)
guidata(hObject, handles);
m_sequence_1=num2str(m_sequence_1);
set(handles.ssequence,'String',m_sequence_1);
guidata(hObject, handles);

% --- Executes on button press in dodajss.
function dodajss_Callback(hObject, eventdata, handles)
% hObject    handle to dodajss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a =get(handles.syg,'String');
guidata(hObject, handles);
b =a-'0';
l = length(a);
k=1;
gold=get(handles.ssequence,'String');
guidata(hObject, handles);
c=gold-'0'
c=c'
g=length(gold);
for i=1:l
for j=1:g
    spread(1,k)=xor(b(1,i),c(1,j));
    k=k+1;
end
end
guidata(hObject, handles);
%spreads=xor(b,c);
hold on
axes(handles.axes3)
l=length(spread)
t=0:1:(l-1);
stem(t,spread);
spread
spread=spread';
spread=num2str(spread);
set(handles.dsequence,'String',spread);
guidata(hObject, handles);

% --- Executes on button press in sekwencja2.
function sekwencja2_Callback(hObject, eventdata, handles)
% hObject    handle to sekwencja2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
register1=[1 1 1 1 1];
code1=zeros(1,31);
for i=1:31
   temp = mod(register1(2)+register1(5),2);
   code1(i) = 2*register1(5)-1;
   for j=5:-1:2
      register1(j)=register1(j-1);
   end
   register1(1) = temp;
end
%m-sequence (5432)
%register2=[1 1 1 1 1];%initial fill
%code2=zeros(1,31);
%for i=1:31
%   temp = mod(register2(2)+register2(3)+register2(4)+register2(5),2);
%   code2(i) = 2*register2(5)-1;
%   for j=5:-1:2
%      register2(j)=register2(j-1);
%   end
%   register2(1) = temp;
%end
m_sequence_1=code1;
%m_sequence_2=code2;
m_sequence_1=m_sequence_1';
%m_sequence_2=m_sequence_2';
m_sequence_1 =m_sequence_1>0;
%m_sequence_2 =m_sequence_2>0;
axes(handles.axes4)
stem(m_sequence_1)
guidata(hObject, handles);
m_sequence_1=num2str(m_sequence_1);
set(handles.osequence,'String',m_sequence_1);
guidata(hObject, handles);

% --- Executes on button press in odejmijawgn.
function odejmijawgn_Callback(hObject, eventdata, handles)
% hObject    handle to odejmijawgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% b=get(handles.dawgn,'String');
% guidata(hObject, handles);
% r=str2num(b);
% r=r';
% c=get(handles.sawgn,'String');
% guidata(hObject, handles);
% n=str2num(c);
% n=n';
% E_s = 1; %Moc sygna³u
% t=0.01:0.01:0.2;
% ds=get(handles.dsequence,'String');
% ds=str2num(ds);
% ds=ds';
% guidata(hObject, handles);
% N=length(ds);
% N
% % index = 1:length(t);
% w=2*pi*5;
% sigma1 = E_s*cos(w*t); %Funkcje bazowe
% sigma2 = E_s*sin(w*t);
% % sigma1
% % sigma2
% % for ii=1:N-1 %Detekcja sygna³u
% % 	x(ii) = sum(sigma1.*r(ii*length(t) + index));
% % 	y(ii) = sum(sigma2.*r(ii*length(t) + index));
% % end
% % x
% % y
% % guidata(hObject, handles);
% % for ii=1:N-1 %Decyzja dla otrzymanych bitów
% %     if(x(ii) < 0)
% %         received_sig(ii) = 1;
% %     else
% %         received_sig(ii) = 2;
% %     end
% % end
% % guidata(hObject, handles);
% m_sig=r-n
% axes(handles.axes8)
% % if and(received_sig(ii)-signal(ii),1)
% % 	plot(X(ii),Y(ii),'ro')
% % else
% % 	plot(X(ii),Y(ii),'bx')
% % end
% %Tworzenie wektorów sygna³u
% for ii=1:2
% 	s(ii,:) = E_s*cos(w*t-2*pi*ii/2);
% end
% %Modulacja losowego strumienia bitów
% signal = m_sig;
% % for i=1:length(signal)
% %     if(m_sig(i) <= 0)
% %         signal(i)=0
% %     else
% %         signal(i)=1
% %     end
% % end
% for i=1:length(m_sig)
%     if(m_sig(i)<0)
%         m_sig(i)=0;
%     else
%         m_sig(i)=1;
%     end
% end
% sig=m_sig+1;
% m_sig = [];
% for ii=1:N-1
% 	m_sig = [m_sig s(sig(ii),:)];
% end
k=get(handles.dawgn,'String');
guidata(hObject, handles);
k=str2num(k);
k=k';
t1=0:0.01:.99;
r1=sin(2*pi*t1);
r2=fliplr(r1);
l=length(k)+length(r2)-1;
d1=fft(k,l);
d2=fft(r2,l);
d=d1.*d2;
p=ifft(d,l);
axes(handles.axes8)
plot(p(1:1500)) %Transmitowany sygna³
guidata(hObject, handles);
set(handles.oawgn,'String',p)
guidata(hObject, handles);

% --- Executes on button press in sygnal.
function sygnal_Callback(hObject, eventdata, handles)
% hObject    handle to sygnal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
spread =get(handles.dmod,'String');
guidata(hObject, handles);
spread=str2num(spread);
los = length(spread);
spread=spread';
m=get(handles.syg,'String');
m=str2num(m);
m=m';
c=get(handles.osequence,'String');
c=str2num(c);
c=c';
i=1;
k=1;
while(k<los)
s=0;
for j=1:length(c)
    temp(1,j)=xor(spread(1,k),c(1,j));
    k=k+1;
    s=s+temp(1,j);
end
if(s==0)
    b2(1,i)=0;
else
    b2(1,i)=1;
end
i=i+1;
end
despreaded_signal=b2;
t = linspace(0,length(spread),length(spread)*100);
t1=0:0.01:.99;
y=int(vpa(despreaded_signal),'t1',0,1);
%b=rectpulse(a(1,:),l)
hold on
axes(handles.axes6)
stem(despreaded_signal)
despreaded_signal=num2str(despreaded_signal);
set(handles.dane,'String',despreaded_signal);


function syg_Callback(hObject, eventdata, handles)
% hObject    handle to syg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of syg as text
%        str2double(get(hObject,'String')) returns contents of syg as a double


% --- Executes during object creation, after setting all properties.
function syg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to syg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ssequence_Callback(hObject, eventdata, handles)
% hObject    handle to ssequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ssequence as text
%        str2double(get(hObject,'String')) returns contents of ssequence as a double


% --- Executes during object creation, after setting all properties.
function ssequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ssequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in demodulator.
function demodulator_Callback(hObject, eventdata, handles)
% hObject    handle to demodulator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% r=get(handles.oawgn,'String');
% r=str2num(r);
% r=r';
% E_s = 1; %Moc sygna³u
% w=2*pi*5;
% t=0.01:0.01:0.2;
% SNR = 0.5;
% E_s = 1; %Moc sygna³u
% No = E_s^2/SNR; %Moc szumu
% sig_n = sqrt(No/2); %Odchylenie standardowe szumu
% w=2*pi*5;
% t=0.01:0.01:0.2;
% ds=get(handles.dsequence,'String')
% ds=str2num(ds);
% guidata(hObject, handles);
% signal=r
% signal=signal'
% for i=1:length(signal)
%     if(signal(i)<0)
%         signal(i)=0;
%     else
%         signal(i)=1;
%     end
% end
% signal=signal+1;
% N=length(ds);
% for ii=1:2
% 	s(ii,:) = E_s*cos(w*t-2*pi*ii/2);
% end
% s=s';
% m_sig = [];
% for ii=1:N-1
% 	m_sig = [m_sig s(signal(ii),:)];
% end
% axes(handles.axes9)
% % for i=1:length(m_sig)
% %     if(m_sig(i)<0)
% %         m_sig(i)=0;
% %     else
% %         m_sig(i)=1;
% %     end
% % end
p=get(handles.oawgn,'String');
p=str2num(p);
p=p';
x=get(handles.dsequence,'String');
guidata(hObject, handles);
x=str2num(x);
% x=x';
% for j=length(x)
%     q(j)=p(100*j);
%     if(q(j)>0)
%         n(j)=1;
%     else
%         n(j)=0;
%     end
% end
% N=length(x);
% t=0.01:0.01:N;
% c=sin(2*pi*t);
% y=p.*c;
dane=get(handles.dsequence,'String');
dane=str2num(dane);
SNR = 50;                                     %10*log10(Eb_N0_dB)+10*log10(1/1); % multiple Eb/N0 values; log10(bitrate/Bandwidth)
NRZ = 2*dane-1;                                     % Kodowanie NRZ
T = 1;                                              % Czas trwania bitu (okres)
fc = 3/T;                                           % Czestotlowosc nosna
t = linspace(0,length(dane),length(dane)*100);      % probki, czasy
N = length(t);                                      % Liczba probek
nosna = sin(pi*fc*t);
bpsk_mod=get(handles.dawgn,'String');
bpsk_mod=str2num(bpsk_mod);
bpsk_mod=bpsk_mod';
for i=1:length(dane)
    odfiltr(i)=bpsk_mod(i).*nosna(i);
end
% for i=1:length(odfiltr)
%     if odfiltr(i)>0
%         bpsk_demod(i)=1;
%     else
%         bpsk_demod(i)=0;
%     end
% end
% szum_p=get(handles.dawgn,'String');
% szum_p=str2num(szum_p);
% kwadrat2 = szum_p.^2;                               % podnosze sygnal do kwadratu 
% filtr=fir1(30,0.99);                                % tworzenie filtru
% pofiltr=filter(filtr,1,kwadrat2);                   % filtrowanie sygnalu 
% x = pofiltr;
% odfiltr=num2str(odfiltr);
% for j=1:length(odfiltr)
%     s=odfiltr(j);
%     y(j)=int(sym('s'),'t',0,1);
% end
% for i=1:length(y)
%     if(y(i)<=0)
%         y(i)=0;
%     else
%         y(i)=1;
%     end
% end
for i=1:length(odfiltr)
    if(odfiltr(i)<0)
        odfiltr(i)=0;
    else
        odfiltr(i)=1;
    end
end
y=odfiltr;
y(1)=0;
axes(handles.axes9);
% plot(t,odfiltr) %Transmitowany sygna³
stem(y);
odfiltr=num2str(y);
guidata(hObject, handles);
set(handles.dmod,'String',y);
guidata(hObject, handles);

% --- Executes on button press in modulator.
function modulator_Callback(hObject, eventdata, handles)
% hObject    handle to modulator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% SNR = 0.5;
% E_s = 1; %Moc sygna³u
% No = E_s^2/SNR; %Moc szumu
% sig_n = sqrt(No/2); %Odchylenie standardowe szumu
% w=2*pi*5;
% t=0.01:0.01:0.2;
% ds=get(handles.dsequence,'String')
% guidata(hObject, handles);
% signal=str2num(ds);
% signal=signal'
% N=length(ds) %Liczba bitów
% % for ii=1:2
% % 	s(ii,:) = E_s*cos(w*t-2*pi*ii/2);
% % end
% % m_sig = [];
% % for ii=1:N
% % 	m_sig = [m_sig s(signal(ii),:)];
% % end
% for i=1:length(signal)
%     if(signal(i)==0)
%         s(i)=sqrt(2*E_s)*cos(w*t);
%     else
%         s(i)=sqrt(2*E_s)*(-cos(w*t));
%     end
% end
% x=get(handles.dsequence,'String');
% guidata(hObject, handles);
% x=str2num(x);
% x=x';
% N=length(x);
% x(x==0)=-1;
% t=0.01:0.01:N;
% c=2*sin(2*pi*t);
% for i=1:1:N
%     m((i-1)*100+1:i*100)=x(i);
% end
% y=c.*m;
dane=get(handles.dsequence,'String');
dane=str2num(dane);
SNR = 50;                                     %10*log10(Eb_N0_dB)+10*log10(1/1); % multiple Eb/N0 values; log10(bitrate/Bandwidth)
NRZ = 2*dane-1;                                     % Kodowanie NRZ
T = 1;                                              % Czas trwania bitu (okres)
fc = 3/T;                                           % Czestotlowosc nosna
t = linspace(0,length(dane),length(dane)*100);      % probki, czasy
N = length(t);                                      % Liczba probek
Lpnb = N/length(dane);                              % Liczba probek na bit
dod_dane = repmat(dane',1,Lpnb);                    % powtarzamy bity Lpnb razy
dod_dane2 = repmat(NRZ',1,Lpnb);
d_d = dod_dane';                                    % Transpozycja wierszy i kolumn
d_d2 = dod_dane2';
d_d = d_d(:)';                                      % Konwersja do jednego wiersza (100x1, 100x0, itd..)
d_d2 = d_d2(:)';     
nosna = sin(pi*fc*t);                               % Fala nosna 
bpsk_mod = d_d2.*nosna;                             % Modulowany przebieg
y=bpsk_mod;
axes(handles.axes7)
plot(t(1:1500),y(1:1500)) %Transmitowany sygna³
guidata(hObject, handles);
y=num2str(y);
set(handles.mod,'String',y);
guidata(hObject, handles);

% --- Executes on button press in dodajawgn.
function dodajawgn_Callback(hObject, eventdata, handles)
% hObject    handle to dodajawgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% b=get(handles.mod,'String');
% guidata(hObject, handles);
% m_sig=str2num(b);
% m_sig=m_sig';
% SNR = 0.5;
% n = awgn(m_sig,SNR); %Wektor szumu
% r = m_sig + n; %Otrzymany sygna³ = s + n
% axes(handles.axes5)
% N=124
% t=1:N-1
y=get(handles.mod,'String');
guidata(hObject, handles);
y=str2num(y);
y=y';
fs=100;
Ej=1;
snr=50;
r=sqrt(fs*Ej*10^(-snr/10))*randn(1,length(y));
%r=randn(1,length(y));
r=r';
k=y+r;
x=get(handles.dsequence,'String');
guidata(hObject, handles);
x=str2num(x);
x=x'
N=length(x);
t=0.01:0.01:N;
axes(handles.axes5)
plot(t(1:1500),k(1:1500)) %Transmitowany sygna³
guidata(hObject, handles);
set(handles.dawgn,'String',k);
guidata(hObject, handles);
set(handles.sawgn,'String',r);
guidata(hObject, handles);

function dsequence_Callback(hObject, eventdata, handles)
% hObject    handle to dsequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dsequence as text
%        str2double(get(hObject,'String')) returns contents of dsequence as a double


% --- Executes during object creation, after setting all properties.
function dsequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dsequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mod_Callback(hObject, eventdata, handles)
% hObject    handle to mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mod as text
%        str2double(get(hObject,'String')) returns contents of mod as a double


% --- Executes during object creation, after setting all properties.
function mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dawgn_Callback(hObject, eventdata, handles)
% hObject    handle to dawgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dawgn as text
%        str2double(get(hObject,'String')) returns contents of dawgn as a double


% --- Executes during object creation, after setting all properties.
function dawgn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dawgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function oawgn_Callback(hObject, eventdata, handles)
% hObject    handle to oawgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of oawgn as text
%        str2double(get(hObject,'String')) returns contents of oawgn as a double


% --- Executes during object creation, after setting all properties.
function oawgn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to oawgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dmod_Callback(hObject, eventdata, handles)
% hObject    handle to dmod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dmod as text
%        str2double(get(hObject,'String')) returns contents of dmod as a double


% --- Executes during object creation, after setting all properties.
function dmod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dmod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function osequence_Callback(hObject, eventdata, handles)
% hObject    handle to osequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of osequence as text
%        str2double(get(hObject,'String')) returns contents of osequence as a double


% --- Executes during object creation, after setting all properties.
function osequence_CreateFcn(hObject, eventdata, handles)
% hObject    handle to osequence (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sawgn_Callback(hObject, eventdata, handles)
% hObject    handle to sawgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sawgn as text
%        str2double(get(hObject,'String')) returns contents of sawgn as a double


% --- Executes during object creation, after setting all properties.
function sawgn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sawgn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ber.
function ber_Callback(hObject, eventdata, handles)
% hObject    handle to ber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
y=get(handles.dane,'String');
guidata(hObject, handles);
y=str2num(y);
y=y';
dane=get(handles.syg,'String');
guidata(hObject, handles);
dane=str2num(dane);
dane=dane';
w=real(y)>0;
z=zeros(1,length(dane));
j=z;
bit=z;
for a=0:length(dane)-1
    for b=1:15
        if(w(100*a+6*b)==0)
            z(a+1)=z(a+1)+1'
        else
            j(a+1)=j(a+1)+1;
        end
    end
end
for a=0:length(dane)-1
    for b=1:15
        if(z(a+1)>j(a+1))
            bit(a+1)=0;
        else
            bit(a+1)=1;
        end
    end
end
nErr=size(find([dane-bit]),2);
nErr
simBER=nErr/length(dane);
simBER
figure();
SNR = -4:12;
EbN0=SNR(1):0.001:SNR(length(SNR));
error = (1/2)*erfc(sqrt(10.^(EbN0/10)));
%sError=
semilogy(50,simBER,'*r',EbN0 ,error,'green');
title('BER')
xlabel('SNR [dB]');
ylabel('BER')



function dane_Callback(hObject, eventdata, handles)
% hObject    handle to dane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dane as text
%        str2double(get(hObject,'String')) returns contents of dane as a double


% --- Executes during object creation, after setting all properties.
function dane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
