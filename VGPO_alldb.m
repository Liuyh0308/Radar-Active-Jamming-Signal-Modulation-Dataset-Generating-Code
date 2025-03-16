%% 速度拖引干扰VGPO

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
folder_path_img = 'D:\\Radar_Jamming_Signal_Dataset\\Trainning_data\\dataset_img\\All_dB\\VGPO';

% 创建文件夹（如果不存在）
if ~exist(folder_path_img, 'dir')
    mkdir(folder_path_img);
end

% 定义序列保存的文件夹路径
folder_path_seq = 'D:\\Radar_Jamming_Signal_Dataset\\Trainning_data\\dataset_seq\\All_dB\\VGPO';

% 创建文件夹（如果不存在）
if ~exist(folder_path_seq, 'dir')
    mkdir(folder_path_seq);
end

index = 1;

for SNR = -20:2:10
    for a = 1:data_num
        %% 目标回波（加入随机延时），并引入假目标，速度拖引VGPO
    
        As = 10^(SNR/20);  % 目标回波幅度
        Aj = 10^((SNR + JSR)/20);  % 干扰回波幅度
    
        sp = zeros(1, samp_num);  % 初始化没有噪声的信号基底,没有噪声的情况下，sp全为0
        
        range_tar = 1 + round(rand(1, 1) * (samp_num - Nfast - 1));  % 随机偏移量，表示目标信号在噪声信号中的起始位置，并确保range_tar不超过允许的范围 
        
        sp(1 + range_tar : Nfast + range_tar) = sp(1 + range_tar : Nfast + range_tar) + As * lfm(1:Nfast);  % 噪声+目标回波，range_tar相当于就是随机时延偏移
    
    
        % 速度拖引干扰：在频域上对目标信号进行频移
        % 根据拖引速度和加速度计算频移
        % 定义拖引速度和加速度
        v_drag = 100+round(rand(1,1)*200);  % 初始拖引速度（米/秒） 100~300
        a_drag = 50+round(rand(1,1)*450);  % 拖引加速度（米/秒^2） 50~500
        c = 3e8;  % 光速（米/秒）
    
        pri_num = 0+round(rand(1,1)*10);
        v_drag_current = v_drag + a_drag * pri_num * PRI;  % 当前时刻的速度
        fd = 2 * v_drag_current / c * fc * 1e6;  % 当前速度引起的频移（多普勒频移）
    
        vgpo_interference = As * lfm(1:Nfast) .* exp(1j * 2 * pi * fd * (0:Nfast-1) * Ts);
           
    
        sp(1 + range_tar : Nfast + range_tar) = sp(1 + range_tar : Nfast + range_tar) + Aj * vgpo_interference;
    
        J = sp;  % 未加入噪声的VGPO
    
    
    
    
    
        %% 将VGPO信号加入噪声
        sp1 = randn([1, samp_num]) + 1j * randn([1, samp_num]);  % 初始化噪声信号基底
        sp1 = sp1 / std(sp1);  % 标准化噪声信号
        J = sp1 + J;
    
        J = J / max(J);  % 归一化
        J_abs = abs(J);
    
        %% 加入多径衰落
    %     num_paths = 3;  % 多径数量
    %     path_gains = [1, 0.7, 0.4];  % 多径增益
    %     path_delays = [0, 1e-7, 2e-7];  % 多径延迟
    %     J_multipath = zeros(size(J));
    %     for i = 1:num_paths
    %         J_multipath = J_multipath + path_gains(i) * circshift(J, round(path_delays(i) * fs));
    %     end
    %     
    % 
    %    %% 加入多普勒频移
    %     doppler_shift = 100;  % 多普勒频移（赫兹）
    %     t = (0:length(J_multipath)-1) / fs;
    %     J_doppler = J_multipath .* exp(1j * 2 * pi * doppler_shift * t);
    %     J=J_doppler;
    
    
        %% 画图
    %     t_plot = linspace(0,PRI,PRI*fs);
    %     figure(1);
    %     plot(t_plot,real(J));
    %     xlabel('时间 (s)');
    %     ylabel('幅度');
    %     title('VGPO时域');
    % 
    %     figure(2),
    %     f_plot = linspace(-fs/2,fs/2,length(t_plot));
    %     plot(f_plot,fftshift(abs(fft(J))))
    %     xlabel('频率')
    %     title('VGPO频域')
    
        figure(3);
        h = hamming(128);
        [S,F,T]=spectrogram(J,h,127,128,fs);
        f = linspace(-fs/2,fs/2,128);
        imagesc(T,f,abs(fftshift(S,1)));
    %     xlabel('时间')
    %     ylabel('频率')
    %     title('VGPO时频图')
    
    %    
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
