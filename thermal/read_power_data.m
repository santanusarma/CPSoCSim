% READ_POWER_DATA reads dynamic power (at max. speed), active leakage and 
% standby leakage over a range of temperatures for each functional unit of 
% a single core. It then normalizes each power component w.r.t. to the
% maximum value and returns them.
%
%   [P_d, P_a, P_s, T] = READ_POWER_DATA(filename, n_units) reads
%   power data from the file in the string filename, assuming the file has
%   power data for n_units. The normalized power matrices are returned 
%   along with the temperature range.
%   
%   [P_d, P_a, P_s, T] = READ_POWER_DATA(filename, n_units, ignore_LDT) 
%   returns leakage values only for the maximum temperature if ignore_LDT
%   is non-zero.
%
% This is the file format for the power data file. It consists of length(T)
% lines, with each line having the following format
%   DieTemperatureInC\tPDynAtMaxSpeed[Unit1]\tPActiveLeak[Unit1]
%   \tPStandbyLeak[Unit1]\tPDynAtMaxSpeed[Unit2]\tPActiveLeak[Unit2]
%   \tPStandbyLeak[Unit2] ... \tPDynAtMaxSpeed[UnitN]\tPActiveLeak[UnitN]
%   \tPStandbyLeak[UnitN]\n
% where N = n_units.
%
% Ravishankar Rao, Arizona State University
% Created, Nov 13 2007

function [P_d, P_a, P_s, T] = read_power_data(filename, n_units, ...
    ignore_LDT, directory)
 
if (nargin > 3) && directory
    File = dir(filename);
    n_skip = 2; % Omit the entries, . and .. and (if present) .svn
    if strcmp(File(n_skip + 1).name, '.svn')
        n_skip = n_skip + 1;
    end
    n_file = length(File) - n_skip;
    
    P_d = zeros(n_units, n_file);
    P_a = zeros(n_units, 13, n_file);
    P_s = zeros(n_units, 13, n_file);
    if nargin > 2 && ignore_LDT
        P_a = zeros(n_units, n_file);
        P_s = zeros(n_units, n_file);
    end
    
    Dl = load('power.txt');
    
    j = 1;
    for i = (n_skip+1):length(File)
        file = ['./power_data/' File(i).name];
        
        % Get raw data
        D = load(file);
        
        % Extract temperature range
        if (nargin > 2) && ignore_LDT
            T = D(end, 1);
        else
            T = Dl(:, 1); % Column 1 of D
        end
        
        % Extract dynamic power data
        C_d = 1 + 1:3:(3*n_units); % Columns 2, 5, 8, 11, ...
        % Same data repeats in all row, so, only pick the first
        P_d(:,j) = D(1, C_d)';
        
        % Extract static and inactive power data, which are interleaved
        C_a = 1 + C_d; % Columns 3, 6, 9, ...
        C_s = 1 + C_a; % Columns 4, 7, 12, ...
        
        if (nargin > 2) && ignore_LDT
            % Only pick leakage values corresponding to the highest temp.
            P_a(:, j) = D(end, C_a)';
            P_s(:, j) = D(end, C_s)';
        else
            % P_a(i, j) = static power of unit i at temp. T(j)
            P_a(:, :, j) = Dl(:, C_a)';
            % P_s(i, j) = inactive power of unit i at temp. T(j)
            P_s(:, :, j) = Dl(:, C_s)';
        end
        j = j + 1;
    end
    % Normalize power data
    P_d = P_d/max(sum(P_d));
    if (nargin > 2) && ignore_LDT
        P_a = P_a/max(sum(P_a));
        P_s = P_s/max(sum(P_s));
    else
        P_a_max = squeeze(P_a(:,end,:));
        P_a = P_a/max(sum(P_a_max));
        P_s_max = squeeze(P_s(:,end,:));
        P_s = P_s/max(sum(P_s_max));
    end
else
    %% Get raw data
    D = load(filename);
    
    %% Extract temperature range
    if (nargin > 2) && ignore_LDT
        T = D(end, 1);
    else
        T = D(:, 1); % Column 1 of D
    end
    
    %% Extract dynamic power data
    C_d = 1 + 1:3:(3*n_units); % Columns 2, 5, 8, 11, ...
    P_d = D(1, C_d)'; % Same data repeats in all row, so, only pick the first
    
    %% Extract static and inactive power data, which are interleaved
    C_a = 1 + C_d; % Columns 3, 6, 9, ...
    C_s = 1 + C_a; % Columns 4, 7, 12, ...
    
    if (nargin > 2) && ignore_LDT
        % Pick leakage values corresponding to the highest temperature only
        P_a = D(end, C_a)';
        P_s = D(end, C_s)';
    else
        P_a = D(:, C_a)'; % P_s(i, j) = static power of unit i at temp. T(j)
        P_s = D(:, C_s)'; % P_i(i, j) = inactive power of unit i at temp. T(j)
    end
    %% Normalize power data
    P_d = P_d/sum(P_d); % sum(P_d) now equals 1
    % Each column sum is less than 1, except the alast column sum, which = 1
    P_a = P_a/sum(P_a(:,end));
    P_s = P_s/sum(P_s(:,end));
end

