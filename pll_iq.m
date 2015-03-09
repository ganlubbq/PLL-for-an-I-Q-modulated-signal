clc
clear all
close all
%Loop Parameters%
digital_bandwidth= 2*pi/100;
damping_factor=sqrt(2)/2;

ki= (4*digital_bandwidth*digital_bandwidth)/...
    (1+2*damping_factor*digital_bandwidth+digital_bandwidth*digital_bandwidth);
kp= (4*damping_factor*digital_bandwidth)/...
    (1+2*damping_factor*digital_bandwidth+digital_bandwidth*digital_bandwidth);

%Input Accumulator%
normalized_freq_1=0.05; normalized_freq_2=0.1;
phase_series=[normalized_freq_1*ones(1,100) normalized_freq_2*ones(1,100)];

%Phase Detector Phase Error%
ph=filter(1,[1 -1],phase_series);
input_complex_sinusoid= exp(1i*2*pi*ph);

ph_sv = zeros(1,200);
ph_sv_1= zeros(1,200);
phase_series_sv = zeros(1,200);

int_reg = 0;
accumulator = 0;
 
%Loop Filter and Phase Accumulator%
for n=1:200
         product =input_complex_sinusoid(n)*exp(-1i*2*pi*accumulator);
         ph_1 = angle(product)/(2*pi);
         ph_sv(n) = ph_1;
         
        int_reg = int_reg + ki*ph_1;
        normalized_freq = kp*ph_1+int_reg;
        phase_series_sv(n) = normalized_freq;
 
        ph_sv_1(n) = accumulator;
        accumulator = accumulator+normalized_freq;
end
output_complex_sinusoid=exp(1i*2*pi*ph_sv_1);
%Phase Error%
phase_error=ph-ph_sv_1;
figure(1)
subplot(3,1,1)
plot(phase_error)  
grid on
title('Phase Detector Phase Error, Input to Loop Filter')
xlabel('Time Index, Output Clock')
ylabel('Phase Error')
 

subplot(3,1,2)
plot(ph_sv)
hold on
plot(phase_series_sv,'r')
hold off
grid on 
title('Loop Filter Output, Input To Phase Accumulator')
xlabel('Time Index, Output Clock')
ylabel('Loop Control')


subplot(3,1,3)
plot(ph)
hold on
plot(ph_sv_1,'r')
hold off
grid on
title('Phase Series of the Two Phase Accumulators')
xlabel('Time Index, Input Clock')
ylabel('Phase Profile')

%Real and Imaginary Parts of i/p and o/p Sinusoids%
s=exp(1i*2*pi*ph_sv_1);
figure(2)
subplot(2,1,1)
plot(real(input_complex_sinusoid))
hold on
plot(real(output_complex_sinusoid),'r')
hold off
grid on
title('Real Part of the i/p and o/p Sinusoids')
xlabel('Time Index')
ylabel('Amplitude')
 
subplot(2,1,2)
plot(imag(input_complex_sinusoid))
hold on
plot(imag(output_complex_sinusoid),'r')
hold off
grid on
title('Imaginary Part of the i/p and o/p Sinusoids')
xlabel('Time Index')
ylabel('Amplitude')

Fs = 4; %Sampling Frequency of the Filter 
Fd = 1; %Sampling Frequency of Digital i/p Signal
delay =10; % Filter Group Delay


r=0.50; 
num=rcosine(Fd,Fs,'sqrt',r,delay); %Nyquist Filter
num_1=num/max(num); %Shaping Filter Scaling 
num_2=num_1/(num_1*num_1');%Proper Scaling for Matched Filter 

%Forming 16-QAM Symbols%
Number_of_Symbols=1000;
Signal_16QAM=(floor(4*rand(1,Number_of_Symbols))-1.5)/1.5+...
    1j*(floor(4*rand(1,Number_of_Symbols))-1.5)/1.5;

%1-to-4 upsampling filters (Shaping Filters)%
h1=reshape([num_1 zeros(1,3)],4,21);
register=zeros(1,20);

x=zeros(1,4*Number_of_Symbols); % Shaping Filter o/p Array
m=0; % o/p Clock Index
for n=1:Number_of_Symbols
    register=[Signal_16QAM(n) register(1:20)]; % Shift i/p Sample into Shaping Filter
    for k=1:4
        x(m+k)=register*h1(k,:)';
    end
    m=m+4;
end

%Eye Diagram%
figure(3)
subplot(2, 1, 1)
plot(0,0)
hold on
for n=(1:4:4*Number_of_Symbols-4)
    plot((-1:1/2:1),real(x(n:n+4)))
end
hold off
grid on
title('Eye Diagram of the Modulated Signal')

%Constellation Diagram%
subplot(2,1,2)
plot(x(1:4:4*Number_of_Symbols),'rx')
grid on
axis('square')
xlabel('Inphase')
ylabel('Quadrature')
title('Constellation Diagram of the Modulated Signal')

sigma=0.02 ;
complex_noise=(sigma*randn(1,4* Number_of_Symbols))+...
    1j*(sigma*randn(1,4* Number_of_Symbols));
 complex_input_sequence=x+complex_noise;  
rotated_sequence= complex_input_sequence.*exp(1j*2*pi/36); %Rotating Complex i/p Sequence
 
%Eye Diagram%
figure(4)
subplot(2, 1, 1)
plot(0,0)
hold on
for n=(1:4:4*Number_of_Symbols-4)
    plot((-1:1/2:1),real(rotated_sequence(n:n+4)))
end
hold off
grid on
title('Eye Diagram of the Rotated Signal')

%Constellation Diagram%
subplot(2,1,2)
plot(rotated_sequence(1:4:4*Number_of_Symbols),'rx')
grid on
axis('square')
xlabel('Inphase')
ylabel('Quadrature')
title('Constellation Diagram of the Rotated Signal')

y=conv(rotated_sequence,num_2);

%Eye Diagram%
figure(5)
subplot(2, 1, 1)
plot(0,0)
hold on
for n=(1:4:4*Number_of_Symbols-4)
    plot((-1:1/2:1),real(y(n:n+4)))
end
hold off
grid on
title('Eye Diagram of the Matched Filter Output Signal')

%Constellation Diagram%
subplot(2,1,2)
plot(y(1:4:4*Number_of_Symbols),'rx')
grid on
axis('square')
xlabel('Inphase')
ylabel('Quadrature')
title('Constellation Diagram of the Matched Filter Output Signal')

flag=1;
mf_register=zeros(1,80);

for n=1:4*Number_of_Symbols
slicer_rotated_sequence=rotated_sequence(n);
mf_register=[slicer_rotated_sequence mf_register(1:80)];
slicer_rotated_sequence_1(n)=mf_register*num_2';


if flag==1 %detect count
        slicer_rotated_sequence_1_real=real(slicer_rotated_sequence_1(n));
          slicer_rotated_sequence_1_real_det=1;
          if slicer_rotated_sequence_1_real<0
              slicer_rotated_sequence_1_real_det=-1;
          end
          if abs(slicer_rotated_sequence_1_real)<2/3;
              slicer_rotated_sequence_1_real_det=slicer_rotated_sequence_1_real_det/3;
          end
         slicer_rotated_sequence_1_imag=imag(slicer_rotated_sequence_1(n));
          slicer_rotated_sequence_1_imag_det=1;
          if slicer_rotated_sequence_1_imag<0
              slicer_rotated_sequence_1_imag_det=-1;
          end
          if abs(slicer_rotated_sequence_1_imag)<2/3;
              slicer_rotated_sequence_1_imag_det=slicer_rotated_sequence_1_imag_det/3;
          end 
          
          
          sum(n)=slicer_rotated_sequence_1_real_det+1j*slicer_rotated_sequence_1_imag_det;
          
          conjugate_product(n)=slicer_rotated_sequence_1(n)*...
              (slicer_rotated_sequence_1_real_det-1j*slicer_rotated_sequence_1_imag_det);

%ATAN%
output_angle(n)=angle(conjugate_product(n))/(2*pi);
end
flag=flag+1;
   
       if flag==5
         flag=1;   
       end
end
figure(6)
subplot(3,1,1)
plot(slicer_rotated_sequence_1(1:4:end),'x')
axis('square')
hold on 
grid on
plot(sum(1:4:end),'rx')
hold off
axis('square')
xlabel('Inphase')
ylabel('Quadrature')
title('2D Slicer i/p and o/p Constellation');

subplot(3,1,2)
plot(conjugate_product(1:4:end),'x')
grid on
axis('square')
xlabel('Inphase')
ylabel('Quadrature')
title('Conjugate Product Constellation');

subplot(3,1,3) 
plot(output_angle);
grid on 
xlabel('Time')
ylabel('Amplitude')
title('Output Time Series of ATAN')

%Loop Parameters%
digital_bandwidth= 2*pi/1000;
damping_factor=sqrt(2)/2;
%damping_factor=1*damping_factor;
ki= (4*digital_bandwidth*digital_bandwidth)/...
    (1+2*damping_factor*digital_bandwidth+digital_bandwidth*digital_bandwidth);
kp= (4*damping_factor*digital_bandwidth)/...
    (1+2*damping_factor*digital_bandwidth+digital_bandwidth*digital_bandwidth);
phase_accumulator=0;
flag=1;
register_1=zeros(1,80);
int_reg=0;
reg_hold=0;

for n= 1:4*Number_of_Symbols
    rotated_sequence_1=rotated_sequence(n)*exp(-1j*2*pi*phase_accumulator);
    register_1=[rotated_sequence_1 register_1(1:80)];
    rotated_sequence_2(n)=register_1*num_2';
    if flag==1 %detect count
        rotated_sequence_2_real=real(rotated_sequence_2(n));
        rotated_sequence_2_real_det=1;
        if rotated_sequence_2_real<0
            rotated_sequence_2_real_det=-1;
        end
        if abs(rotated_sequence_2_real)<2/3;
            rotated_sequence_2_real_det=rotated_sequence_2_real_det/3;
        end
        rotated_sequence_2_imag=imag(rotated_sequence_2(n));
        rotated_sequence_2_imag_det=1;
        if rotated_sequence_2_imag<0
            rotated_sequence_2_imag_det=-1;
        end
        if abs(rotated_sequence_2_imag)<2/3;
            rotated_sequence_2_imag_det=rotated_sequence_2_imag_det/3;
        end
        
        
        sum(n)=rotated_sequence_2_real_det+1j*rotated_sequence_2_imag_det;
        
        conjugate_product(n)=rotated_sequence_2(n)*...
            (rotated_sequence_2_real_det-1j*rotated_sequence_2_imag_det);
        %ATAN%
        output_angle=angle(conjugate_product(n))/(2*pi);
        int_reg=int_reg+ki*output_angle;
        reg_hold=int_reg+kp*output_angle;
        output_angle_1(n)=output_angle;
    end
    
    filter_reg(n)=reg_hold;
    flag=flag+1;
    if flag==5
        flag=1;
    end
    phase_accumulator=phase_accumulator+reg_hold;
end

figure(7)
subplot(3,1,1)
plot(rotated_sequence_2(1:4:end),'x')
axis('square')
hold on
grid on
plot(sum(1:4:end),'rx')
hold off
axis('square')
xlabel('Inphase')
ylabel('Quadrature')
title('De-Spinned Constellation');

subplot(3,1,2)
plot(conjugate_product(1:4:end),'x')
grid on
axis('square')
xlabel('Inphase')
ylabel('Quadrature')
title('Conjugate Product Constellation');

subplot(3,1,3)
plot(output_angle_1);
grid on
xlabel('Time')
ylabel('Amplitude')
title('Output Time Series of ATAN')

%pause(15);
%close all;
