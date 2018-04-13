% Verify Extracted Kinetics
% Tim Dorn
% 9th Apr 2008
% 
% --------------------------------------------------------------------
% Usage: verifyKineticsVicon(nfigure, S, X, FP, saveImages, markersDyn*, camUp*)
% --------------------------------------------------------------------
% 
% Inputs:   nfigure = number of figures to plot
%           S = GRF with time vector as first row
%           X = CoP with time vector as first row
%           FP = all force plate data (from getKinetics)
%           saveImages = toggles image save (0 = no, 1 = yes)
%           markersDyn* = dynamic marker coordinates structure (from getMarkers)
%               (if none, then markersDyn = [])
%           camUp* = camera up vector (default is +Z direction [0 0 1])
% 
% 
% Notes
% -----
% 
%   Verifies that the extracted data is consistent by plotting the GRF
%   at the CoP on the force plates
% 
%   IMPORTANT: make sure that all input data is in the SAME coordinate
%              system!!!!
%  
% --------------------------------------------------------------------

function verifyKinetics(nfigure, S, X, FP, saveImages, markersDyn, camUp)

usage = 'Usage: verifyKinetics(nfigure, S, X, FP, saveImages, markersDyn*, camUp*)';

if nargin == 5,
    markersDyn = [];
    camUp = [0 0 1];
elseif nargin == 6,
    camUp = [0 0 1];
elseif nargin ~= 7,
    disp(usage)
    return
end

warning off MATLAB:divideByZero
FPcorners = FP.corners;
ntime = length(S(1,:));
scale = markersDyn.divide_to_meters;
e = 0.2;
f = round([1:ntime/(nfigure-1):ntime, ntime]);
loadLabels;

if isstruct(markersDyn),
   r = markersDyn.a2vRatio;
   markerpos = markersDyn.data(:, 2:end);
   l = length(markerpos(1,:));
end


for i=1:6
   co.S(i) = pchip(S(1,:), S(i+1,:));
   co.X(i) = pchip(X(1,:), X(i+1,:));
end


for i = 1:nfigure
   fnum = figure('Name', sprintf('%s: #%d / #%d', strtok(markersDyn.c3dFile, '.'), i, nfigure));

   for j = 1:6
      Si(j) = ppval(co.S(j),f(i));
      Xi(j) = ppval(co.X(j),f(i));
   end
   
   
   for j = 2:3
       
       % plot the GRF at the CoP
       subplot(1, 2, j-1, 'CameraUpVector', camUp)
       quiver3(Xi(1), Xi(2), Xi(3), Si(1)/scale, Si(2)/scale, -Si(3)/scale,...
		        0, 'b-', 'LineWidth', 4); hold on;
       quiver3(Xi(4), Xi(5), Xi(6), Si(4)/scale, Si(5)/scale, -Si(6)/scale,...
		        0, 'r--', 'LineWidth', 4); hold on;
        
        % plot glocal origin point
        plot3(0,0,0, 'ko', 'MarkerSize', 8)
            
        % set some plot properties     
        if camUp(2) == 1,
            if j == 2,
                % xz view
                set(gca, 'CameraUpVector', camUp)
                view(0,0)
            elseif j == 3,
                camPos = [-11.9 9 -8.4];
                set(gca, 'CameraPosition', camPos)
                set(gca, 'CameraUpVector', camUp)
            end
            
        elseif camUp(3) == 1,
            if j == 2,
                % xy view
                view(2)
            elseif j == 3,
                view(3)
            end
        end
        
        axis equal
        grid on
        hold on
        title(sprintf('Frame # %d ,rGRF(Y) = %.2f, lGRF(Y) = %.2f', ...
            f(i), Si(1+2), Si(1+5)));
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
   
                    
        % plot the legend     
        if j == 3, 
            legend('Right Foot', 'Left Foot', 'Location', 'SouthOutside'); 
        end
        
        % plot the force plates
        % Corner Reference: FPcorners(CornerNumber, Coord, PlateNumber)
        cnrEndOrder = [2,3,4,1];
        [x,y,z] = size(FPcorners);
        for plate = 1:z
            for corner = 1:x
                cnrStart = corner;
                cnrEnd = cnrEndOrder(corner);
           
                plot3([FPcorners(cnrStart, 1, plate), FPcorners(cnrEnd, 1, plate)], ...
                    [FPcorners(cnrStart, 2, plate), FPcorners(cnrEnd, 2, plate)], ...
                    [FPcorners(cnrStart, 3, plate), FPcorners(cnrEnd, 3, plate)], ...
                    'Color', [1, 0.6, 1], 'LineWidth', 2)
            end
            text(FPcorners(1, 1, plate), FPcorners(1, 2, plate), ...
                FPcorners(1, 3, plate), ['FP', num2str(FP.order(plate))]);
        end

        % set the plot limits
        if camUp(2) == 1,       % camUp = [0 1 0]
            xlim([min(reshape(squeeze(FPcorners(:,1,:)), x*z, 1)) - e, ...
                max(reshape(squeeze(FPcorners(:,1,:)), x*z, 1)) + e])
            zlim([min(reshape(squeeze(FPcorners(:,3,:)), x*z, 1)) - e, ...
                max(reshape(squeeze(FPcorners(:,3,:)), x*z, 1)) + e])
%             zlim(-[max(reshape(squeeze(FPcorners(:,3,:)), x*z, 1)) + e, ...
%                 min(reshape(squeeze(FPcorners(:,3,:)), x*z, 1)) - e])
            ylim([0, 2])
            
        else                    % camUp = [0 0 1]
            xlim([min(reshape(squeeze(FPcorners(:,1,:)), x*z, 1)) - e, ...
                max(reshape(squeeze(FPcorners(:,1,:)), x*z, 1)) + e])
            ylim([min(reshape(squeeze(FPcorners(:,2,:)), x*z, 1)) - e, ...
                max(reshape(squeeze(FPcorners(:,2,:)), x*z, 1)) + e])
            zlim([0, 2])
        end
   
        % plot dynamic markers (if given as an optional input)
        if isstruct(markersDyn),
           fv = floor(((f(i) - 1)/r) + 1);      % video frame
           
           % get vectors of x, y, z of marker positions in vicon coords
           % for the particular frame in question
           xM = []; yM = []; zM = [];
           for k = 1:l
               if mod(k,3) == 1, xM = [xM ; markerpos(fv,k)]; end
               if mod(k,3) == 2, yM = [yM ; markerpos(fv,k)]; end
               if mod(k,3) == 0, zM = [zM ; markerpos(fv,k)]; end
           end

           % plot 3d markers
           scale = 1000;
           plot3(xM/scale, yM/scale, zM/scale, ...
               'MarkerFaceColor', [0.6, 0.8, 0.2] ,'Marker', 'o', ...
               'LineStyle', 'none', 'MarkerSize', 5, ...
               'Color', [0.4, 0.4, 0.2]);  
       end
       view3d rot 
        
       pos = get(gcf,'Position');
       F(:,i) = getframe(gcf,[0 0 pos(3) pos(4)]);
       
       if j == 2, F1(:,i) = getframe;
       elseif j == 3, F2(:,i) = getframe;
       end
   end
   
   if saveImages == 1,
        num = sprintf('%0.3d', i);
        saveas(fnum, sprintf('%sverifyKinetics%s_%s.emf', ...
            glab.infoDirectory, num, markersDyn.c3dFile), 'emf');
   end
   figs(i) = gcf;

end

% F1 = movie looking downward
% F2 = movie with 3D isometric view
figure('Name', [markersDyn.c3dFile, ': Movie Frame'], 'NumberTitle', 'off')
playmovie(F, ceil(nfigure/4), 1);
playmovie(F, ceil(nfigure/6), 1, 'slider');
% playmovie(F2, ceil(nfigure/4), 1);
% playmovie(F2, ceil(nfigure/6), 1, 'slider');
close(figs)

warning on MATLAB:divideByZero

