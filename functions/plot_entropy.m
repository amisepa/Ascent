function plot_entropy(entropyData, chanlocs, chanIdx)

chanlabels = {chanlocs(chanIdx).labels};
x = [ chanlocs(chanIdx).X ]';
y = [ chanlocs(chanIdx).Y ]';
z = [ chanlocs(chanIdx).Z ]';

% Rotate X Y Z coordinates
% rotate = 0;       %nosedir = +x
rotate = 3*pi/2;    %nosedir = +y
% rotate = pi;      %nosedir = -x
% rotate = pi/2;
allcoords = (y + x.*sqrt(-1)).*exp(sqrt(-1).*rotate);
x = imag(allcoords);
y = real(allcoords);

% Project 3D positions on 2D plane if not already done
chanpos(:,1) = x;
chanpos(:,2) = y;
chanpos(:,3) = z;

if all(chanpos(:,3)==0)
    coord = chanpos(:,1:2); % positions already projected on a 2D plane
else
    coord = chanpos; % use 3-D data for plotting
end

% deal with near-0 values
idx = entropyData < 0.01;
if sum(idx) > 0
    entropyData(idx) = 0.0001;
    warning(['Channel ' chanlocs(idx).labels ' is probably a bad channel.'])
end

% Scalp topo
if var(entropyData) > 0.1
    figure('color','w');
    topoplot(entropyData, chanlocs, 'emarker', {'.','k',15,1},'electrodes','labels');
    clim([min(entropyData) max(entropyData)]);
    c = colorbar; colormap('parula'); %ylabel(c,'Power spectral difference (log)','FontSize',12)
    % title('Entropy','FontSize',10); 
    c.Label.String = 'Entropy';
    c.Label.FontSize = 12;
    c.Label.FontWeight = 'normal';
else
    disp('Not enough variance across electrodes to plot an informative scalp topography.')
end

% 3D figure allowing to open entropy data for each electrode
p = figure('color','w');
% p.Position = [100 100 540 400];
axis equal
axis vis3d
axis off
hold on

% adj = mean(entropyData(1,:))*5; % to scale marker size

for iChan = 1:size(entropyData,1)
    
    if length(entropyData(iChan,:)) == 1 % measures with one value per channel
        % 3D plot of entropy values at electrode locations
        p(iChan) = plot3(coord(iChan,1),coord(iChan,2),coord(iChan,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor', 'k', ...
            'Marker','o','MarkerSize',5);
        
        % Display channel label + entropy value for each channel
        text(coord(iChan,1)-15,coord(iChan,2)+10,coord(iChan,3), ...
           sprintf('%s: %6.1f',chanlabels{iChan}, ...
           round(entropyData(iChan,:),2)),'FontSize',10,'fontweight','bold');

    else % for multiscales, take area under the curve as sensor size
        p(iChan) = plot3(coord(iChan,1),coord(iChan,2),coord(iChan,3), ...
            'MarkerEdgeColor','k','MarkerFaceColor', 'k', ...
            'Marker','o','MarkerSize', 5, 'UserData',iChan, ...
            'ButtonDownFcn', @(~,~,~) buttonCallback(entropyData(iChan,:), coord(iChan,:), chanlabels{iChan}));

        % Display channel label above each electrode
        text(coord(iChan,1)-7,coord(iChan,2)+10,coord(iChan,3), ...
            sprintf('%s %6.3f',chanlabels{iChan}), ...
            'FontSize',10,'fontweight','bold');
        title('[Click on sensors to display entropy values]', ...
            'Position', [1 120 1], 'fontweight', 'bold')

    end
end

%% subfunction to display entropy values in the plot where use clicks

function buttonCallback(tmpdata, coor, label)

% Entropy measures with only one value per channel
figure('color','w','Position', [500 500 280 210]);
plot(tmpdata,'linewidth',2,'color','black'); % blue: [0, 0.4470, 0.7410]
% area(tmpdata,'linewidth',2);
title(label,'FontSize',14)
% xticks(2:nScales); xticklabels(join(string(scales(:,2:end)),1)); xtickangle(45)
% xlim([2 nScales]);
xlabel('Time scale','FontSize',12,'fontweight','bold'); 
ylabel('Entropy','FontSize',12,'fontweight','bold')

% end
