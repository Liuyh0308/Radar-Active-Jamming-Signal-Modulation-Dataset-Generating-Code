%% 间歇采样转发干扰 ISRJ

clear
clc
close all

%% 常数
data_num=200;   %干扰样本数

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

JSR=5;   %干信比dB

%% 生成LFM发射信号
Nr  = round(fs*PRI);  %一帧信号的长度,即在一个脉冲重复周期PRI内，采样得到的信号样本总数
Nfast = round(fs*taup); % 计算脉宽长度，即在一个脉冲宽度（taup）内，采样得到的LFM信号样本总数，也就是LFM有效存在的时间（采样个数来计）
t = (-Nfast/2:(Nfast/2 - 1))*Ts;  % 脉内的时间序列
lfm = [exp(1j*(2*pi*fc*t+pi*k*(t).^2)), zeros(1, Nr-Nfast)];  %LFM信号 复包络 并将其扩展到一帧信号长度Nr的长度
       
samp_num = Nr; %距离窗点数

% 定义图像和序列的根文件夹路径
root_folder_img = 'D:\\Radar_Jamming_Signal_Dataset\\Test_data\\dataset_img';
root_folder_seq = 'D:\\Radar_Jamming_Signal_Dataset\\Test_data\\dataset_seq';


% 循环生成不同SNR的信号数据集
for SNR = -20:2:10
    % 定义图像保存的文件夹路径
    folder_path_img = sprintf('%s\\%d_dB\\ISRJ', root_folder_img, SNR);

    % 创建图像保存文件夹（如果不存在）
    if ~exist(folder_path_img, 'dir')
        mkdir(folder_path_img);
    end

    % 定义序列保存的文件夹路径
    folder_path_seq = sprintf('%s\\%d_dB\\ISRJ', root_folder_seq, SNR);

    % 创建序列保存文件夹（如果不存在）
    if ~exist(folder_path_seq, 'dir')
        mkdir(folder_path_seq);
    end


    for a=1:data_num
       %%  目标回波（加入随机延时）
    
        As=10^(SNR/20);%目标回波幅度
        Aj=10^((SNR+JSR)/20);%干扰回波幅度
    
        sp = zeros(1, samp_num);  %初始化没有噪声的信号基底,没有噪声的情况下，sp全为0
        range_tar = 1 + round(rand(1, 1) * (samp_num - Nfast-2000 - 1)); % 随机偏移量，表示目标信号在噪声信号中的起始位置，并确保range_tar不超过允许的范围 
        sp(1 + range_tar : Nfast + range_tar) = sp(1 + range_tar : Nfast + range_tar) + As * lfm(1:Nfast); % 噪声+目标回波，range_tar相当于就是随机时延偏移
        echo=sp;  %未加入噪声的回波信号
        
        %% ISRJ，采样+转发
        repetion_times=[4,2,1];   %重复转发次数
        period=[4e-6,5e-6,10e-6];    %表示采样次数，taup / period =5, 4或2次，采样脉冲周期
        duty=[20,25,33.33,50];  %占空比
     
        index1=1+round(rand(1,1)* (length(period) - 1));  %生成随机索引，用于从 period, duty 和 repetion_times 数组中选择参数
        index2=1+round(rand(1,1)* (length(repetion_times) - 1));
        
        period1=period(index1);
        duty1=duty(index2);
        repetion_times1=repetion_times(index2);
    
    
        m =round(taup / period1) ; %采样次数
        n =repetion_times1 ; %转发次数
        len = Nfast/m/n;%每次采样的长度,m*n需为可被Nfast=2000整除的整数
        len_m = len*n ; %第n次转发拒采样时的长度
    
        delay_time=period1*duty1*0.01;  %延迟一个采样脉冲时间，即采样后立即转发
        delay_num=ceil(delay_time / (1/fs));  %ceil()为进一法取整，表示一个延迟时间内程序采样点数
        
        len_t = delay_num;%不同采样转发间隔点数
    
        J = zeros(1,Nr); %此变量用于存储间歇采样转发干扰ISRJ信号
        
        for i=1:m
            %采样
            temp = echo((1+range_tar+(i-1)*len_m) : (1+range_tar+len+(i-1)*len_m-1));
            for j=1:n
                %转发
                J(1+range_tar+(i-1)*len_m + (j-1)*len+(i-1)*len_t : 1+range_tar +j*len+(i-1)*len_m-1+(i-1)*len_t) = Aj*temp * 2;  %间歇采样转发干扰  ,PS.. (* 1.5是为了让ISRJ信号更明显一点)
            end
        end
       %% 将ISRJ信号加入噪声
        sp1=randn([1,samp_num])+1j*randn([1,samp_num]);  %初始化噪声信号基底
        sp1=sp1/std(sp1); %标准化噪声信号
        sp1(1+range_tar:1+range_tar +n*len+(m-1)*len_m-1+(m-1)*len_t)=...
        sp1(1+range_tar:1+range_tar +n*len+(m-1)*len_m-1+(m-1)*len_t)+J(1+range_tar:1+range_tar +n*len+(m-1)*len_m-1+(m-1)*len_t); %加入噪声的ISRJ干扰信号
        J=sp1;  %J为最终的加入噪声的回波ISRJ
    
%        %% 加入多径衰落
%         num_paths = 3;  % 多径数量
%         path_gains = [1, 0.7, 0.4];  % 多径增益
%         path_delays = [0, 1e-7, 2e-7];  % 多径延迟
%         J_multipath = zeros(size(J));
%         for i = 1:num_paths
%             J_multipath = J_multipath + path_gains(i) * circshift(J, round(path_delays(i) * fs));
%         end
%         
%     
%        %% 加入多普勒频移
%         doppler_shift = 100;  % 多普勒频移（赫兹）
%         t1 = (0:length(J_multipath)-1) / fs;
%         J_doppler = J_multipath .* exp(1j * 2 * pi * doppler_shift * t1);
%         J=J_doppler;
      
        %% 画图
        J=J/max(J); %归一化
        J_abs=abs(J);
    
    
    %     t_plot = linspace(0,PRI,PRI*fs);
    %     figure(1);
    %     plot(t_plot,real(J));
    %     xlabel('时间 (s)');
    %     ylabel('幅度');
    %     title('ISRJ时域');
    % 
    %     figure(2),
    %     f_plot = linspace(-fs/2,fs/2,length(t_plot));
    %     plot(f_plot,fftshift(abs(fft(J))))
    %     xlabel('频率')
    %     title('ISRJ频域')
    
        figure(3);
        h = hamming(128);
        [S,F,T]=spectrogram(J,h,127,128,fs);
        f = linspace(-fs/2,fs/2,128);
        imagesc(T,f,abs(fftshift(S,1)));
    %     xlabel('时间')
    %     ylabel('频率')
    %     title('ISRJ时频图')
    
    
    
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
        J_fft = fftshift(fft(J(1 + range_tar : Nfast + range_tar), 1024));       % 计算长度为1024的FFT
        file_name_seq = fullfile(folder_path_seq, sprintf('%d.mat', a));
        save(file_name_seq, 'J_fft');       % 保存为 .mat 文件
    
    end
end
