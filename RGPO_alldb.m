%% 距离拖引干扰RGPO

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
folder_path_img = 'D:\\Radar_Jamming_Signal_Dataset\\Trainning_data\\dataset_img\\All_dB\\RGPO';

% 创建文件夹（如果不存在）
if ~exist(folder_path_img, 'dir')
    mkdir(folder_path_img);
end

% 定义序列保存的文件夹路径
folder_path_seq = 'D:\\Radar_Jamming_Signal_Dataset\\Trainning_data\\dataset_seq\\All_dB\\RGPO';

% 创建文件夹（如果不存在）
if ~exist(folder_path_seq, 'dir')
    mkdir(folder_path_seq);
end

index = 1;

for SNR = -20:2:10
    for a = 1:data_num
        %% 目标回波（加入随机延时），并引入假目标，距离拖引RGPO
        As = 10^(SNR/20);      % 目标回波幅度
        Aj = 10^((SNR+JSR)/20); % 干扰回波幅度
    
        sp = zeros(1, samp_num);  % 初始化没有噪声的信号基底,没有噪声的情况下，sp全为0
        
        range_tar = 1 + round(rand(1, 1) * (samp_num - Nfast - 1-4000)); % 随机偏移量，表示目标信号在噪声信号中的起始位置，并确保range_tar不超过允许的范围 
        
        sp(1 + range_tar : Nfast + range_tar) = sp(1 + range_tar : Nfast + range_tar) + As * lfm(1:Nfast); % 噪声+目标回波，range_tar相当于就是随机时延偏移
        
        % 初始化干扰信号位置
        % 初始化距离波门的位置
        range_gate_start = 500;   % 距离波门起始位置
        range_gate_width = 50;   % 距离波门宽度（样本数）
        pull_off_rate = 100+round(rand(1,1)*400);       % 每个PRI周期后移的距离（样本数）  拖引速度（样本/秒）：  pull_off_rate/PRI = 100~400m/s
    
        range_jamming = range_gate_start;
    
        K = 1+round(rand(1,1)*5);
    
        for pri = 1:K   % 相当于设置了拖引时间 K个PRI时间
            
            range_jamming = range_jamming + pull_off_rate;
        end
        sp(1 + range_tar + range_jamming : Nfast + range_tar + range_jamming) = sp(1 + range_tar + range_jamming : Nfast + range_tar + range_jamming) + Aj * lfm(1:Nfast); 
        J = sp;  % 未加入噪声的RGPO
     
        %% 将RGPO信号加入噪声
        sp1 = randn([1, samp_num]) + 1j * randn([1, samp_num]);  % 初始化噪声信号基底
        sp1 = sp1 / std(sp1); % 标准化噪声信号
        J = sp1 + J;
    
        J = J / max(J); % 归一化
        J_abs = abs(J);
    
        %% 画图
    %     t_plot = linspace(0,PRI,PRI*fs);
    %     figure(1);
    %     plot(t_plot,real(J));
    %     xlabel('时间 (s)');
    %     ylabel('幅度');
    %     title('RGPO时域');
    % 
    %     figure(2),
    %     f_plot = linspace(-fs/2,fs/2,length(t_plot));
    %     plot(f_plot,fftshift(abs(fft(J))))
    %     xlabel('频率')
    %     title('RGPO频域')
    
        figure(3);
        h = hamming(128);
        [S, F, T] = spectrogram(J, h, 127, 128, fs);
        f = linspace(-fs/2, fs/2, 128);
        imagesc(T, f, abs(fftshift(S, 1)));
    %     xlabel('时间')
    %     ylabel('频率')
    %     title('RGPO时频图')
    
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
        J_fft = fftshift(fft(J(1 + range_tar : Nfast + range_tar), 1024));  % 计算长度为1024的FFT
        file_name_seq = fullfile(folder_path_seq, sprintf('%d.mat', index));
        save(file_name_seq, 'J_fft');  % 保存为 .mat 文件

        index = index + 1;
    end

end
