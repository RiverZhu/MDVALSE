% Demo for the MDVALSE algorithm for line spectral signals
% For details, refer to the paper
% Q. Zhang, J. Zhu*, N. Zhang and Z. Xu, Multidimensional Variational Line Spectra Estimation,
% IEEE Signal Processing Letters, vol. 27, 2020 
% This code is written by Qi Zhang and Jiang Zhu. If you have any problems, please feel free to contact
% jiangzhu16@zju.edu.cn.
clc;
clear variables;
rng(123)
D    = 4;   % the number of dimensions;
N    = 5;  
M    = N*ones(1,D); % we set M_1 = M_2 = ... = M_D = N;
K    = 3;   % the number of targets;       
d    = pi/N;        
SNR  = 10;          

%% Generate the line spectral signal
a_theta = @(index,omega)exp((1i*(index-1))*omega);
[c_M{1:D}] = ndgrid(1:M);
D_set_vec = reshape(cat(D,c_M{:}),[],D);

omega = zeros(K,D);
for i = 1:D
    omega(1,i) = pi * (2*rand - 1);
    for k = 2:K
        th = pi * (2*rand - 1);
        while min(abs((wrapToPi(th-omega(1:k-1,i))))) < d
            th = pi * (2*rand - 1);
        end
        omega(k,i) = th;
    end
end
a = zeros(prod(M),1);
x = zeros(prod(M),1);
R   = 1 + .2*randn(K,1);                       
w  = R.*exp(1i*2*pi*rand(K,1));                  
for k = 1:K
    a = a_theta(D_set_vec,(omega(k,:))');
    x = x + a*w(k);
end

Pn  = sum(abs(x).^2)/prod(M)*10^(-SNR/10);       
eps = sqrt(0.5*Pn)*(randn(prod(M),1)+1i*randn(prod(M),1));  
y   = x + eps; 

%% Perform MDVALSE algorithm
outinfVALSE = MDVALSE(y, N, x, D, M);
omega_est = outinfVALSE.freqs;
weight_est = outinfVALSE.amps;

%% We display the estimates
mse_vec = outinfVALSE.mse;
mse = mse_vec(end);
disp(['The true frequencies are']);
disp(num2str(round(sort(omega),2)))
disp(['The estimated frequencies are']);
disp(num2str(round(sort(omega_est),2)))
disp(['The mse of the line spectral signal is ', num2str(round(10*log10(mse),2)), ' dB.']);

%% If D == 2, we plot the figure

if D == 2
    figure(1)
    set(gca, 'FontSize', 16);             
    stem3(omega(:,1),omega(:,2),abs(w), 'ro-', 'Linewidth',2,'MarkerSize',10)
    hold on
    stem3(omega_est(:,1),omega_est(:,2),abs(outinfVALSE.amps), 'kx--', 'Linewidth',2,'MarkerSize',10)
    xlabel('Dimension 1', 'FontSize', 16);
    ylabel('Dimension 2', 'FontSize', 16);
    zlabel('Amplitude', 'FontSize', 16);
    xlim([-pi, pi])
    ylim([-pi, pi])
    legend('Ground Truth','MDVALSE')
end





