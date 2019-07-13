% Define figures for real time update

% Start figure handle to view fields in time and space
f1 = figure; %figure('Position', [1000, 100, 500, 100]);
% ax(3) = subplot(3,1,3);
power = zeros(ydim,xdim);
% h = imagesc(delta*1*(1:1:xdim),(delta*1*(1:1:ydim))',power,[-80,0]);  
% h = imagesc(delta*1*(1:1:xdim),(delta*1*(1:1:ydim))',power,[0 0.001]); 
% h = imagesc(delta*1*(1:1:xdim),(delta*1*(1:1:ydim))',power,[-1e-0 1e-0]);
% h = imagesc(delta*1*(1:1:xdim),(delta*1*(1:1:ydim))',power,[-1e-3 1e-3]);
% h = imagesc(delta*1*(1:1:xdim),(delta*1*(1:1:ydim))',power,[-1e-5 1e-5]);
h = imagesc(delta*1*(1:1:xdim),(delta*1*(1:1:ydim))',power,[-1e-1 1e-1]);
set(h.Parent,'FontSize',15)
% h = imagesc(delta*1*(1:1:xdim),(delta*1*(1:1:ydim))',power); 
xlabel('Horizontal direction (m)','FontSize',20); ylabel('Depth direction (m)','FontSize',20); cbarlabel = colorbar; ylabel(cbarlabel,'Electric Field','FontSize',20);%ylabel(cbarlabel,'Power (dB)')
% title(['2D TM_z, with sampling = \lambda / ' num2str(cellsperwavelength)]);
axis equal
axis([0 x_dist 0 y_dist])
set(f1,'Color',[1 1 1])
% linkaxes(ax);