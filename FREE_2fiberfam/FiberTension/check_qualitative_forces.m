% check_qualitative_forces.m
%   Generate a plot accross all combinations of fiber angles of the 
%   directions of the net force and torque on a FREE at its relaxed
%   configuration.
%
%   Calls solveFB.m to determine values of parameters at test pressure P.

% choose the test pressure
P_test = 1;
% choosing resting parameters for FREE
r_rest = 3/16;
L_rest = 5;

qual = [0 0 0 0];

for i = -44:44
    for j = -44:44
        
        if i == 0
            gama0 = 1;
        else
            gama0 = 2*i;
        end
        
        if j == 0
            betta0 = 1;
        else
            betta0 = 2*j;
        end
        
        
        x_rest = [0, deg2rad(gama0), deg2rad(betta0), r_rest, L_rest, 0];

        [dL, dphi] = qual_FB(P_test, x_rest);
        
        qual = [qual; gama0, betta0, dL, dphi];
    end
end


%% Visualize Results in an accessible way

% Filter out cases where gama is not bigger in magnitude
qplot = zeros(181,3);
count = 1;
for k = 1:length(qual(:,1))

    count = count + 1;
    if (abs(qual(k,1)) > abs(qual(k,2)))
        qplot(count,:) = [qual(k,1), qual(k,2), 2*qual(k,3) + 1*qual(k,4)];
    end

%    qplot(qual(k,1) + 91,:) = [qual(k,1), qual(k,2), 200*qual(k,3) + 100*qual(k,4)];

end

% plot the results
M = zeros(181,181);

for m = 1:length(qplot(:,1))

   M(qplot(m,1) + 91, qplot(m,2) + 91) = qplot(m,3);
   
end

figure
[X,Y] = meshgrid(-90:90,-90:90);

% create a colormap having RGB values of dark green,
%light green, white, dark red and light red.
map2 = [1 0 0; 0 1 0; 1 1 1; 1 0.64 0; 0 0 1];
%use the user defined colormap for figure.
colormap(map2);
%plot the figure
shading('faceted');
h = pcolor(X,Y,M');
set(h, 'EdgeColor', 'none');
colorbar('Ticks',[-3,-1,0,1,3],...
         'TickLabels',{'Contrace/Clockwise','Contract/Counter Clockwise','Unknown','Extend/Clockwise','Extend/Counter Clockwise'})


