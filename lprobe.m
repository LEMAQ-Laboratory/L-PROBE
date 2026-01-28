function [spike_vector,Y3D,n,position]= lprobe(X,options)
%
%[spike_vector,Y3D,n,position]= lprobe(X,options);
%
% L-PROBE: Local PROminence-Based Elimination
%
% This algorithm identifies and corrects spikes in Raman imaging data
% using a neighborhood-based local prominence and width filtering strategy.
%
% Developed at LEMAQ - Mass Spectrometry and Chemometrics Laboratory
% Federal University of Paraná (UFPR), Brazil.
%
% -------------------------------------------------------------------------
% INPUTS:
%   X       : Original hyperspectral cube (spiked) after using fsmload.m
%   options : Struct with the following optional fields:
%       .k                    : Number of standard deviations for confidence
%                               limit (default: 5)
%       .prominence_threshold : Threshold for the prominence step (default: 20)
%       .width_threshold      : FWHM width threshold (default: 5)
%       .plots                : Display plots [0 (NONE) or 1 (ALL)]
%
% OUTPUTS:
%   Y3D          : Processed data cube (spike-free)
%   n            : Total number of spikes removed
%   spike_vector : Raman shift positions where spikes were detected
%   position     : Spatial coordinates of spikes in the hyperspectral cube
%
% -------------------------------------------------------------------------
% Developed by: Prof. Dr. Frederico Luis Felipe Soares
% Contact: frederico.soares@ufpr.br
% -------------------------------------------------------------------------

% Check options
if exist('options','var')
    if isfield(options,'k')
        k = options.k;
    else
        k = 5;
    end
    if isfield(options,'prominence_threshold')
        prominence_thrd = options.prominence_threshold;
    else
        prominence_thrd = 20;
    end
    if isfield(options,'width_threshold')
        width_threshold = options.width_threshold;
    else
        width_threshold = 5;
    end
    if isfield(options,'plots')
        plots = options.plots;
    else
        plots = 0;
    end
%else
   % k = 5;
    %prominence_thrd = [];
    %width_threshold = 20;
    %lambda = 10e-2;
    %p = 0.01;
    %plots = 1;
    
end

%initial time
tic

%dimensions recognition
[m1,m2,m3]=size(X);
if m3 == 1
    error('A hypercube is required! It must have X x Y x Wv dimensions')
else
    X3D=X;
    X2D=reshape(X,(m1*m2),m3);
    Y3D = X3D;
end
[m1,m2,m3]=size(X3D);

position = zeros(m1, m2, m3);%
% Spike detection based on FWHM and Proeminence
spike_vector = [];
for i=1:size(X2D,1)
    data = X2D(i,:);
    
    if isempty(prominence_thrd)
        [~,P] = ecdf(data);
        prominence_threshold = 3*std(P(1:round(length(P)*.25)));
    else
        prominence_threshold = prominence_thrd;
    end        
     
    [spikes_pks, spikes_locs, spikes_widths, spikes_prom] = findpeaks(...
        data, ...
        'MinPeakProminence', prominence_threshold, ...
        'MaxPeakWidth', width_threshold);

    % Ratio Prominence/FWHM
    min_ratio = 5; % Spikes have high ratio
    valid_spikes = (spikes_prom ./ spikes_widths) > min_ratio;
    spikes_locs = spikes_locs(valid_spikes);
    spikes_pks = spikes_pks(valid_spikes);
    
    % Create the prominence vector
    mask = ~ismember(spikes_locs, spike_vector);
    spike_vector = [spike_vector, spikes_locs(mask)];
end

%counter spikes
n=0;

f = waitbar(0,'1','Name','Merging pixels...');
% 
% [Xasls,~]=asls_baseline(X2D, lambda, p);% Test the ASLS before the routine 
% X3D = reshape(Xasls, m1, m2, m3);

for wv=1:size(spike_vector,2) %slices
    c = spike_vector(wv);
    %c = proeminence(spec)
    
    % Update waitbar and message
    waitbar(wv/size(spike_vector,2),f,['Slice ' num2str(wv) ' of ' num2str(size(spike_vector,2))])
    
    %--------------------------------------------
    %situation 1 (Inside the image)
    % Define neighborhood indices
    a = 2:(m1-1);
    a1 = a-1;
    a2 = a+1;
    b = 2:(m2-1);
    b1 = b-1;
    b2 = b+1;
    
    % Get all 8-neighborhood values at once
    neighbors = cat(4, ...
        X3D(a1,b1,c), X3D(a,b1,c), X3D(a2,b1,c), ...
        X3D(a1,b,c),                X3D(a2,b,c), ...
        X3D(a1,b2,c), X3D(a,b2,c), X3D(a2,b2,c));
    
    % Sort along the 4th dimension
    sorted_neighbors = sort(neighbors, 4);
    
    % Compute statistics
    ym = mean(sorted_neighbors(:,:,:,2:7), 4);
    yd = std(sorted_neighbors(:,:,:,2:7), 0, 4);
    
    % Create masks
    mask = X3D(a,b,c) < (ym + k*yd);
    
%    
%     % Create original neighbors (Must think a better way)
%     neighbors_original = cat(4, ...
%         Y3D(a1,b1,c), Y3D(a,b1,c), Y3D(a2,b1,c), ...
%         Y3D(a1,b,c),               Y3D(a2,b,c), ...
%         Y3D(a1,b2,c), Y3D(a,b2,c), Y3D(a2,b2,c));
%     
%     % Sort along the 4th dimension
%     sorted_neighbors_original = sort(neighbors_original, 4);
%     
    % Update Y3D and position
    Y3D(a,b,c) = mask .* X3D(a,b,c) + ~mask .* mean(sorted_neighbors(:,:,:,2:6), 4);
%     Y3D(a,b,c) = mask .* Y3D(a,b,c) + ~mask .* mean(sorted_neighbors_original(:,:,:,2:6), 4);
    position(a,b,c) = ~mask;
    
    % Count replacements
    n = n + sum(~mask(:));
    
    %--------------------------------------------
    %situation 2.1 (vertical edges of image)
    a = 2:(m1-1);
    a1 = a-1;
    a2 = a+1;
    for b_edge = [1, m2]
        
        if b_edge == 1
            b = 1;
            b2 = b+1;  % Left edge
            
            % Get all 5-neighborhood values at once
            neighbors = cat(4, ...
                X3D(a1,b,c), X3D(a1,b2,c), ...
                X3D(a,b2,c),             ...
                X3D(a2,b,c), X3D(a2,b2,c));
        else
            b = m2;
            b1 = b-1; % Right edge
            
            % Get all 5-neighborhood values at once
            neighbors = cat(4, ...
                X3D(a1,b1,c), X3D(a1,b,c), ...
                              X3D(a, b1,c), ...
                X3D(a2,b1,c), X3D(a2,b,c));
        end
        
        % Sort along the 4th dimension
        sorted_neighbors = sort(neighbors, 4);
        
        % Compute statistics
        ym = mean(sorted_neighbors(:,:,:,2:4), 4);
        yd = std(sorted_neighbors(:,:,:,2:4), 0, 4);
        
        % Create masks
        mask = X3D(a,b,c) < (ym + k*yd);
        
        % Update Y3D and position
        Y3D(a,b,c) = mask .* X3D(a,b,c) + ~mask .* mean(sorted_neighbors(:,:,:,2:4), 4);
        position(a,b,c) = ~mask;
        
        % Count replacements
        n = n + sum(~mask(:));
    end
    
    %--------------------------------------------
    %situation 2.2 (horizontal edges of image)
    b = 2:(m2-1);
    b1 = b-1;
    b2 = b+1;
    
    for a_edge = [1, m1]
        
        if a_edge == 1
            a = 1;
            a2 = a+1;  % up edge
            
            % Get all 5-neighborhood values at once
            neighbors = cat(4, ...
                X3D(a, b1,c),              X3D(a2,b2,c), ...
                X3D(a2,b1,c), X3D(a2,b,c), X3D(a2,b2,c));
            
        else
            a = m1;
            a1 = a-1;  % down edge
            
            % Get all 5-neighborhood values at once
            neighbors = cat(4, ...
                X3D(a1,b1,c), X3D(a1,b,c), X3D(a1,b2,c), ...
                X3D(a, b1,c),              X3D(a ,b2,c));
        end
        
        % Sort along the 4th dimension
        sorted_neighbors = sort(neighbors, 4);
        
        % Compute statistics
        ym = mean(sorted_neighbors(:,:,:,2:4), 4);
        yd = std(sorted_neighbors(:,:,:,2:4), 0, 4);
        
        % Create masks
        mask = X3D(a,b,c) < (ym + k*yd);
        
        % Update Y3D and position
        Y3D(a,b,c) = mask .* X3D(a,b,c) + ~mask .* mean(sorted_neighbors(:,:,:,2:4), 4);
        position(a,b,c) = ~mask;
        
        % Count replacements
        n = n + sum(~mask(:));
        
    end
    
    
    
    %--------------------------------------------
    %situation 3 (Corners of the image)
    
    for a_edge = [1, m1]
        for b_edge = [1, m2]            
            if a_edge == 1
                if b_edge == 1 % up right corner
                    a=1;
                    b=1;
                    a2=a+1;
                    b2=b+1;
                    
                    % Get all 3-neighborhood values at once
                    neighbors = cat(4, ...
                                      X3D(a2,b ,c), ...
                        X3D(a ,b2,c), X3D(a2,b2,c));
                    
                else  % up left corner
                    a=1;
                    b=m2;
                    a2=a+1;
                    b1=b-1;
                    
                    % Get all 3-neighborhood values at once
                    neighbors = cat(4, ...
                        X3D(a2,b1,c),             ...
                        X3D(a2,b1,c), X3D(a2,b ,c));
                end
            else
                if b_edge == 1 % down right corner
                    a=m1;
                    b=1;
                    a1=a-1;
                    b2=b+1;
                    
                    % Get all 3-neighborhood values at once
                    neighbors = cat(4, ...
                        X3D(a1,b ,c), X3D(a1,b2,c), ...
                                      X3D(a ,b2,c));
                    
                else  % down left corner
                    a=m1;
                    b=m2;
                    a1=a-1;
                    b1=b-1;
                    
                    % Get all 3-neighborhood values at once
                    neighbors = cat(4, ...
                        X3D(a1,b1,c), X3D(a1,b ,c),...
                        X3D(a ,b1,c)              );
                end
            end
            
            % Sort along the 4th dimension
            sorted_neighbors = sort(neighbors, 4);
            
            % Compute statistics
            ym = mean(sorted_neighbors(:,:,:,:), 4);
            yd = std(sorted_neighbors(:,:,:,:), 0, 4);
            
            % Create masks
            mask = X3D(a,b,c) < (ym + k*yd);
            
            % Update Y3D and position
            Y3D(a,b,c) = mask .* Y3D(a,b,c) + ~mask .* mean(sorted_neighbors(:,:,:,:), 4);
            position(a,b,c) = ~mask;
            
            % Count replacements
            n = n + sum(~mask(:));
            
        end
    end
end
close(f)


% Sort locs from spike
spike_vector = sort(spike_vector);

%%PLOTS
if options.plots == 1
    %plot spikes free
    Y2D=reshape(Y3D,(m1*m2),m3);
    figure;
    subplot (2,1,1);plot(X2D');
    subplot (2,1,2);plot(Y2D');
    figure; plot(Y2D');
    
    %plot of spike identification in the matrix (data)
    figure;
    plot(X2D', 'b-'); hold on; plot(spike_vector, max(X2D(:,spike_vector)),'or');

end


disp(['Number of replaced spikes ' num2str(n)]);
%final time
elapsedTime = toc;
disp(['Elapsed time of analysis ' num2str(elapsedTime)]);