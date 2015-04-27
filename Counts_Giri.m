function counts

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "The preRC protein ORCA/LRWD1 Organizes Heterochromatin by 
% Assembling Histone H3 Lysine 9 Methyltransferases on Chromatin"
%
% Sumanprava Giri, Vasudha Aggarwal, Julien Pontis, Zhen Shen , 
% Arindam Chakraborty, Abid Khan, Craig Mizzen, Kannanganattu V. 
% Prasanth, Slimane Ait-Si-Ali, Prof. Taekjip Ha, Supriya Prasanth
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Count the average number of fluorescence from all the images in a
% given folder. The input files are .traces files obtained after 
% identification of single-molecules using IDL scripts. traces files 
% include the intensity values for every frame for every identified 
% single molecule. IDL scripts are available at bio.physics.illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all;
fclose('all');

%% ===== Input Directory with .traces files ==========================
path = input('Directory with .traces files :   ');
cd(path);
Number_Files = dir;  % all files in the current directory defined by path
B = pwd;
D = findstr(B,'\');
E = max(D);
CurrentDirectory = B(E+1:end);

% =========== Define the matric to store average intensity values =======
avg_intensity_temp = zeros(0,1);
combine_avg_intensity = zeros(0,0);
total_intensity = zeros(0,0);
count_number = zeros(0,0);

%% ============== read only the .traces files =======================
for i = 1:size(Number_Files)
    if Number_Files(i).isdir == 0
        s = Number_Files(i).name;
        if s(end-5:end) == 'traces'
            disp(s);
            file_id = fopen(s,'r');
            length = fread(file_id,1,'int32');    % length of the .traces
            fprintf('Frame length of movie: %d\n', length);
            number_of_molecules = fread(file_id,1,'int16');
            fprintf('Number of single-molecules in the movie: %d\n', number_of_molecules/2);
            raw = fread(file_id, number_of_molecules*length,'int16');
            disp('Done reading data.');
            
            %=======separate 'raw' into the individual single molecule data =====
            index = (1:number_of_molecules*length);
            Data = zeros(number_of_molecules,length);
            intensity = zeros(number_of_molecules/2,length);
            Data(index) = raw(index);
            for i=1:(number_of_molecules/2)
                intensity(i,:) = Data(i*2-1,:); 
                avg_intensity(i) = sum(intensity(i,2:11),2)/10; % avg over 10 frames
                avg_intensity_temp = [avg_intensity_temp avg_intensity(i)];
            end
            combine_avg_intensity = [combine_avg_intensity avg_intensity_temp];
            final_avg_intensity = combine_avg_intensity';
        end
    end  
end    
time = (0:(length-1))*0.1;  % time vector, 0.1 sec is frame rate
            
%% ================= Plot the intensity histogram =========================
figure;
hist(final_avg_intensity(:,1),80);
grid on; zoom on;
title(s);

%% ==================== Determine Intensity Cut offs ====================
fcutoff1 = input('low cutoff intensity: ','s');
cutoff1 = str2num(fcutoff1);
fcutoff2 = input('high cutoff intensity: ','s');
cutoff2 = str2num(fcutoff2);
sel_index = (final_avg_intensity > cutoff1) & (final_avg_intensity < cutoff2);
sel_avg_intensity = final_avg_intensity(sel_index,:);           
save([CurrentDirectory '_Intensity_' num2str(cutoff1) '_' num2str(cutoff2) '.dat'],'sel_avg_intensity','-ascii');
close;

B = pwd;
n_count = 0;

%% ======================== read the data again ======================
fname = [CurrentDirectory '_counts.dat'];
for i = 1:size(Number_Files)
    if Number_Files(i).isdir == 0
        s = Number_Files(i).name;
        if s(end-5:end) == 'traces'
            file_id = fopen(s,'r');
            length = fread(file_id,1,'int32');    % length of the .traces
            number_of_molecules = fread(file_id,1,'int16');
            raw = fread(file_id, number_of_molecules*length,'int16');
            
            %=======separate 'raw' into the individual single molecule data =====
            index = (1:number_of_molecules*length);
            Data = zeros(number_of_molecules,length);
            intensity = zeros(number_of_molecules/2,length);
            Data(index) = raw(index);
            for i=1:(number_of_molecules/2)
                intensity(i,:) = Data(i*2-1,:);
                cutoff_value = mean(intensity(i,(2:5)));
                if cutoff_value > cutoff1 & cutoff_value < cutoff2
                    n_count = n_count + 1;
                end
            end
            fid_temp = fopen(fname, 'a+');
            temp_output=[s ['      ',int2str(n_count)]];
            fprintf(fid_temp,'%s\n',temp_output);
            fclose(fid_temp);
            count_number = [count_number n_count];
            n_count = 0;
        end
    end
end

%% ============== Calculate Avg and Std of the count ==================
Ave_number_Total = mean(count_number);
SD_Total = std(count_number);
fid_temp = fopen(fname, 'a+');
fprintf(fid_temp,'\n\n');
temp_output=['Average Total is:       ' ];
fprintf(fid_temp,'%s',temp_output);
fprintf(fid_temp,'%g\n',Ave_number_Total);
temp_output=['Standard deviation is:   '];
fprintf(fid_temp,'%s',temp_output);
fprintf(fid_temp,'%g\n',SD_Total);
fclose(fid_temp);

Ave_number_round = round(Ave_number_Total);
SD_round = round(SD_Total);

figure;
hist(sel_avg_intensity(:,1),80);
zoom on;
temp=axis;
title ([CurrentDirectory '     Average = ' num2str(Ave_number_round) '\pm' num2str(SD_round)])
text((temp(2))*0.7,(temp(4)-temp(3))*7/8, ['Int. Range = ' num2str(cutoff1) '-' num2str(cutoff2)], 'FontSize',10, 'Color','blue');

end
