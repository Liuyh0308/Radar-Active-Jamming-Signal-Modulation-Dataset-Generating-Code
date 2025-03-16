%% 噪声乘积式灵巧噪声干扰 MNJ
clear
clc
close all

%% 常数
data_num = 800;   % 干扰样本数

%% 信号源参数
fs = 100e6;  % 采样频率 100MHz
Ts = 1/fs;  % 采样间隔

PRI = 100e-6; % 脉冲重复周期 100μs
fc = 10e6;   % 载频

%% LFM信号参数
B = 40e6;  % 信号带宽 40MHz
taup = 40e-6; % 信号脉宽 40μs
k = B/taup; % 调频斜率

JSR = 5;  % 干信比dB

%% 生成LFM发射信号
Nr = round(fs*PRI);  % 一帧信号的长度
Nfast = round(fs*taup); % 计算脉宽长度
t = (-Nfast/2:(Nfast/2 - 1))*Ts;  % 脉内的时间序列
lfm = [exp(1j*(2*pi*fc*t+pi*k*(t).^2)), zeros(1, Nr-Nfast)];  % LFM信号

samp_num = Nr; % 距离窗点数

% 定义图像保存的文件夹路径
folder_path_img = 'D:\\Radar_Jamming_Signal_Dataset\\Trainning_data\\dataset_img\\All_dB\\MNJ';

% 创建文件夹（如果不存在）
if ~exist(folder_path_img, 'dir')
    mkdir(folder_path_img);
end

% 定义序列保存的文件夹路径
folder_path_seq = 'D:\\Radar_Jamming_Signal_Dataset\\Trainning_data\\dataset_seq\\All_dB\\MNJ';

% 创建文件夹（如果不存在）
if ~exist(folder_path_seq, 'dir')
    mkdir(folder_path_seq);
end

index = 1;

for SNR = -20:2:10
    for a = 1:data_num
       %%  目标回波（加入随机延时）
    
        As=10^(SNR/20);%目标回波幅度
        Aj=10^((SNR+JSR)/20);%干扰回波幅度
    
        sp = zeros(1, samp_num);  %初始化没有噪声的信号基底,没有噪声的情况下，sp全为0
        range_tar = 1 + round(rand(1, 1) * (samp_num - Nfast - 1)); % 随机偏移量，表示目标信号在噪声信号中的起始位置，并确保range_tar不超过允许的范围 
        sp(1 + range_tar : Nfast + range_tar) = sp(1 + range_tar : Nfast + range_tar) + As * lfm(1:Nfast); % 噪声+目标回波，range_tar相当于就是随机时延偏移
        echo=sp;  %未加入噪声的回波信号
    
       %% 将NP信号加入噪声
        sp1=randn([1,samp_num])+1j*randn([1,samp_num]);  %初始化噪声信号基底
        sp1=sp1/std(sp1); %标准化噪声信号
        echo_with_noise=sp1+Aj*echo * 2;    % * 2是为了让NP信号更明显一点
    
      %% 噪声乘积干扰是窄带高斯白噪声与回波信号进行相乘
        L=Nr; %噪声样本数量
        wgn_noise=wgn(L,1,0);  % 生成宽带高斯白噪声
    
        %% 设计带通滤波器，让宽带噪声变为窄带白噪声
        bw_filter=[5e6,7e6,9e6,11e6]; %噪声带宽选择
        index1=1+round(rand(1,1)* (length(bw_filter) - 1));
    
        f0 = fc;  % 噪声中心频率=lfm载频
        bw = bw_filter(index1);   % 噪声带宽
        d = fdesign.bandpass('N,F3dB1,F3dB2', 10, f0 - bw / 2, f0 + bw / 2, fs);%指定滤波器设计的参数：滤波器阶数 N=10 和 3dB 截止频率 F3dB1、F3dB2。fs 是采样频率，确保滤波器设计与信号的采样频率一致
        Hd = design(d, 'butter');  % 设计巴特沃斯带通滤波器
    
        %% 通过滤波器生成窄带高斯白噪声
        narrowband_noise = filter(Hd, wgn_noise);
    
        %% 对窄带噪声进行希尔伯特变换
        analytic_narrowband_noise = hilbert(narrowband_noise);
    
        %% 计算噪声乘积干扰信号
        J = analytic_narrowband_noise.'.*echo_with_noise;
        
    
         %% 画图
        J=J/max(J); %归一化
        J_abs=abs(J);
    
    
    %     t_plot = linspace(0,PRI,PRI*fs);
    %     figure(1);
    %     plot(t_plot,real(J));
    %     xlabel('时间 (s)');
    %     ylabel('幅度');
    %     title('NP时域');
    % 
    %     figure(2),
    %     f_plot = linspace(-fs/2,fs/2,length(t_plot));
    %     plot(f_plot,fftshift(abs(fft(J))))
    %     xlabel('频率')
    %     title('NP频域')
    
        figure(3);
        h = hamming(128);
        [S,F,T]=spectrogram(J,h,127,128,fs);
        f = linspace(-fs/2,fs/2,128);
        imagesc(T,f,abs(fftshift(S,1)));
    %     xlabel('时间')
    %     ylabel('频率')
    %     title('NP时频图')
    
    
        axis off;  % 关闭坐标轴
        set(gca, 'Position', [0 0 1 1]);  % 去掉图像周围的空白
    
    
         % 保存图像
        frame = getframe(gcf);
        img = frame.cdata;
        resized_img = imresize(img, [224, 224]);
        file_name_img = fullfile(folder_path_img, sprintf('%d.png', index));
        imwrite(resized_img, file_name_img);
        close;    % 关闭当前图形窗口
    
    
        % 保存序列
        J_fft = fftshift(fft(J(1 + range_tar : Nfast + range_tar), 1024));       % 计算长度为1024的FFT
        file_name_seq = fullfile(folder_path_seq, sprintf('%d.mat', index));
        save(file_name_seq, 'J_fft');       % 保存为 .mat 文件

        index = index + 1;
    end

end