%% 雷达LFM线性调频信号
clear
clc
close all

%% 常数
data_num=10;   %干扰样本数

%% 信号源参数
fs = 100e6;  % 采样频率 100MHz
Ts = 1/fs ;  % 采样间隔


PRI = 100e-6; % 脉冲重复周期 100μs (两个连续脉冲信号之间的时间间隔)
% PRF = 1/PRI;  % 脉冲重复频率
fc  = 10e6;   % 载频

%% LFM信号参数
B = 40e6;  % 信号带宽 40MHz
taup = 40e-6; % 信号脉宽（脉冲宽度） 40μs 一个脉冲的持续时间，即信号从开始到结束的时间长度
k = B/taup; % 调频斜率

JSR=5;  %信噪比dB

%% 生成LFM发射信号
Nr  = round(fs*PRI);  %一帧信号的长度,即在一个脉冲重复周期PRI内，采样得到的信号样本总数
Nfast = round(fs*taup); % 计算脉宽长度，即在一个脉冲宽度（taup）内，采样得到的LFM信号样本总数，也就是LFM有效存在的时间（采样个数来计）
t = (-Nfast/2:(Nfast/2 - 1))*Ts;  % 脉内的时间序列
lfm = [exp(1j*(2*pi*fc*t+pi*k*(t).^2)), zeros(1, Nr-Nfast)];  %LFM信号 复包络 并将其扩展到一帧信号长度Nr的长度

      
samp_num = Nr; %距离窗点数

% 定义图像和序列的根文件夹路径
root_folder_img = 'D:\\Radar_Jamming_Signal_Dataset\\Test_data\\dataset_img';
root_folder_seq = 'D:\\Radar_Jamming_Signal_Dataset\\Test_data\\dataset_seq';

% num_label = 0; %标签设置
% label=zeros(1,data_num)+num_label;  %标签数据,此干扰标签为0

% 循环生成不同SNR的信号数据集
for SNR = -20:2:10
    % 定义图像保存的文件夹路径
    folder_path_img = sprintf('%s\\%d_dB\\LFM', root_folder_img, SNR);

    % 创建图像保存文件夹（如果不存在）
    if ~exist(folder_path_img, 'dir')
        mkdir(folder_path_img);
    end

    % 定义序列保存的文件夹路径
    folder_path_seq = sprintf('%s\\%d_dB\\LFM', root_folder_seq, SNR);

    % 创建序列保存文件夹（如果不存在）
    if ~exist(folder_path_seq, 'dir')
        mkdir(folder_path_seq);
    end
    

    for a=1:data_num
        % 目标回波＋噪声
        sp=randn([1,samp_num])+1j*randn([1,samp_num]);  %噪声信号基底
        sp=sp/std(sp); %标准化噪声信号
        As=10^(SNR/20);%目标回波幅度
        Aj=10^((SNR+JSR)/20);%干扰回波幅度
        range_tar = 1 + round(rand(1, 1) * (samp_num - Nfast - 1)); % 随机偏移量，表示目标信号在噪声信号中的起始位置，并确保range_tar不超过允许的范围 
        
        lfm = [exp(1j*(2*pi*fc*t+pi*k*(t).^2)), zeros(1, Nr-Nfast)];
    
        sp(1 + range_tar : Nfast + range_tar) = sp(1 + range_tar : Nfast + range_tar) + As * lfm(1:Nfast); % 噪声+目标回波，range_tar相当于就是随机时延偏移
        lfm_echo=sp;
        lfm_echo=lfm_echo/max(lfm_echo); %归一化
        lfm_echo_abs=abs(lfm_echo);

%         %% 加入多径衰落
%         num_paths = 3;  % 多径数量
%         path_gains = [1, 0.7, 0.4];  % 多径增益
%         path_delays = [0, 1e-7, 2e-7];  % 多径延迟
%         J_multipath = zeros(size(lfm_echo));
%         for i = 1:num_paths
%             J_multipath = J_multipath + path_gains(i) * circshift(lfm_echo, round(path_delays(i) * fs));
%         end
%         
%     
%        %% 加入多普勒频移
%         doppler_shift = 100;  % 多普勒频移（赫兹）
%         t1 = (0:length(J_multipath)-1) / fs;
%         J_doppler = J_multipath .* exp(1j * 2 * pi * doppler_shift * t1);
%         lfm_echo=J_doppler;
    
        %% 无噪声回波
%         As=10^(SNR/20);%目标回波幅度
%         sp = zeros(1, samp_num);  %初始化没有噪声的信号基底,没有噪声的情况下，sp全为0
%         range_tar = 1 + round(rand(1, 1) * (samp_num - Nfast - 1)); % 随机偏移量，表示目标信号在噪声信号中的起始位置，并确保range_tar不超过允许的范围 
%         sp(1 + range_tar : Nfast + range_tar) = sp(1 + range_tar : Nfast + range_tar) + As * lfm(1:Nfast); % 噪声+目标回波，range_tar相当于就是随机时延偏移
%         lfm_echo=sp;
        
        %% 画图
        t_plot = linspace(0,PRI,PRI*fs);
        figure(1);
        plot(t_plot,real(lfm_echo));
        xlabel('时间 (s)');
        ylabel('幅度');
        title('LFM信号时域');
    
        figure(2),
        f_plot = linspace(-fs/2,fs/2,length(t_plot));
        plot(f_plot,fftshift(abs(fft(lfm_echo))))
        xlabel('频率')
        title('LFM信号频域')
    
        figure(3);
        h = hamming(128);
        [S,F,T]=spectrogram(lfm_echo,h,127,128,fs);
        f = linspace(-fs/2,fs/2,128);
        imagesc(T,f,abs(fftshift(S,1)));
        xlabel('时间')
        ylabel('频率')
        title('LFM信号时频图')
        

        axis off;  % 关闭坐标轴
        set(gca, 'Position', [0 0 1 1]);  % 去掉图像周围的空白
    
    
         % 保存图像
        frame = getframe(gcf);
        img = frame.cdata;
        resized_img = imresize(img, [224, 224]);
        file_name_img = fullfile(folder_path_img, sprintf('%d.png', a));
        imwrite(resized_img, file_name_img);
        close;    % 关闭当前图形窗口
    
    
        % 保存序列
        lfm_echo_fft = fftshift(fft(lfm_echo(1 + range_tar : Nfast + range_tar), 1024));       % 在干扰信号存在区间计算长度为1024的FFT
        file_name_seq = fullfile(folder_path_seq, sprintf('%d.mat', a));
        save(file_name_seq, 'lfm_echo_fft');       % 保存为 .mat 文件
    end

end
