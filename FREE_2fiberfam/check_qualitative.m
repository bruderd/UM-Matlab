% Check qualitative behavior of FREEs. This is an older version that calls
% the control optimizer so it runs slow and generates poor results. For
% better and faster version see check_qualitative_forces.m

Num = 5;
T = 1;
qual = [0 0 0 0];


% for i = -17:17
%     for j = -17:17
for i = 1:10
    for j = 1:10
        L0 = 5;
        gama0 = 8*i + 1;
        betta0 = 8*j + 2;
        x_rest = [0.0001, deg2rad(gama0), deg2rad(betta0), 3/16, L0, 0]';
        [x_star, u_star] = main_2fiberfam_func(x_rest, Num, T);
        
        close all
        
        qual = [qual; gama0, betta0, sign(x_star(5,end)-L0), sign(x_star(6,end))];
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

[X,Y] = meshgrid(-90:90,-90:90);

% create a colormap having RGB values of dark green,
%light green, white, dark red and light red.
map2 = [1 0 0; 0 1 0; 1 1 1; 1 0.64 0; 0 0 1];
%use the user defined colormap for figure.
colormap(map2);
%plot the figure
shading('faceted');
h = pcolor(X,Y,M);
set(h, 'EdgeColor', 'none');
colorbar('Ticks',[-3,-1,0,1,3],...
         'TickLabels',{'Contrace/Clockwise','Contract/Counter Clockwise','Unknown','Extend/Clockwise','Extend/Counter Clockwise'})


