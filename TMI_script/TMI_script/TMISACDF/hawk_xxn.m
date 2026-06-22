function output = hawk_xxn(input, levels, neg_mode, group_by_time)
%   input - 输入图像堆栈 [height×width×frames]
%   levels      - 处理层级数 (正整数)
%   neg_mode    - 负值处理模式: 'abs'或'separate'
%   group_by_time  - 是否按时间排序: true或false

%   output - 处理后的图像堆栈

[height, width, num_frames] = size(input);
output_cells = {};

for l = 0:levels-1
    kernel_width = 2^(l+1);      % 核总宽度
    half_width = 2^l;            % 核半宽
    valid_frames = num_frames - kernel_width + 1;
    
    for s = 1:valid_frames
        % 提取当前核范围帧
        kernel_frames = input(:,:,s:s+kernel_width-1);
        
        % 与小波核卷积
        front_part = sum(kernel_frames(:,:,1:half_width), 3);
        back_part = sum(kernel_frames(:,:,half_width+1:end), 3);
        diff = front_part - back_part;% 卷积值
        
        % 处理负值，俩模式
        if strcmpi(neg_mode, 'abs')
            pos_frame = abs(diff);
            neg_frame = [];
        elseif strcmpi(neg_mode, 'separate')
            pos_frame = max(diff, 0);
            neg_frame = max(-diff, 0);
        end
        
        %结构体 存储每个t时刻，不同尺度下的的小波分解幅值 
        current_frame = struct();
        current_frame.pos = pos_frame;
        current_frame.neg = neg_frame;
        current_frame.center = s + half_width - 1; % 中心帧位置
        
        % 根据排序模式存储
        if group_by_time
            if isempty(output_cells) || length(output_cells) < current_frame.center
                output_cells{current_frame.center} = [];
            end
            output_cells{current_frame.center} = [output_cells{current_frame.center}, current_frame];
        else
            output_cells = [output_cells, {current_frame}];
        end
    end
end

% 整理最终输出堆栈
output = [];
for i = 1:length(output_cells)
    frame_group = output_cells{i};
    for j = 1:length(frame_group)
        output = cat(3, output, frame_group(j).pos);
        if ~isempty(frame_group(j).neg)
            output = cat(3, output, frame_group(j).neg);
        end
    end
end

% 添加微小噪声防止全零
output(1,1,:) = output(1,1,:) + eps;
end