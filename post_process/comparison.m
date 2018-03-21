%clear all
x = h5read('/Volumes/Chest/snapshots_1E7_10_32/snapshots_1E7_10_32_s30.h5','/scales/x/1.0');
z = h5read('/Volumes/Chest/snapshots_1E7_10_32/snapshots_1E7_10_32_s30.h5','/scales/z/1.0');

t = h5read('/Volumes/Chest/snapshots_1E7_10_32/snapshots_1E7_10_32_s30.h5','/scales/sim_time');

dx = x(2)-x(1);


T = h5read('/Volumes/Chest/snapshots_1E7_10_32/snapshots_1E7_10_32_s30.h5','/tasks/b');

T_o = load('/Volumes/Work/Fluids_project/Programs/POD/optimal_solns/optimal_solutions/T_1E7_10.txt');

Nx = 128 ; Ny = 201;

T_o = reshape(T_o,[Ny,Nx]);

y = load('/Volumes/Work/Fluids_project/Programs/POD/optimal_solns/optimal_solutions/y.txt');

alpha = 0.1560021496301813E+002; %1E7, 10

%alpha = 0.4714574564629509E+001; %1E7, 4

Lx = 2*pi/alpha;
x_o = -Lx/2:Lx/(Nx):Lx/2 - Lx/(Nx);

Mz = size(T,1);
Mx = size(T,2);

T_s = T(:,:,1) - repmat(z,[1,Mx]);

y_l = Mz/2;

threshold = 0.16;
    
[T_max,locs] = findpeaks(T_s(y_l,:).*(T_s(y_l,:)>threshold));

%E_l2(1:50,1:size(locs,2))=0;

for j = 1:50
    
    T_s = T(:,:,j) - repmat(z,[1,Mx]);
    
    [T_max,locs] = findpeaks(T_s(y_l,:).*(T_s(y_l,:)>threshold));
    
    % figure(1)
    % for i = 1:size(locs,2)
    %     contourf(x_o + x(locs(i)),-y,T_o,15)
    %     hold on
    % end
    
    % colormap(jet)
    % axis equal tight
    
    %contour(x,z,T_s,15,'LineWidth', 2)
    %colormap(bone)
    
    
    [X,Y] = meshgrid(x_o',-y);
    
    loc = Lx/2;
    
    i_r = find(x < (loc + dx/2) & x > (loc - dx/2));
    
    i_l = find(x < (-loc + dx/2) & x > (-loc - dx/2));
    
    [Xq,Yq] = meshgrid(x(i_l:i_r),z);
    
    x_d = x(i_l:i_r);
    
    T_i = interp2(X,Y,T_o,Xq,Yq,'spline');
    
    figure(1)
    plot(locs(1),'b*')
    hold on
    
    for k = 1:size(locs,2)
        figure(1)
        plot(j,locs(k),'b*')
        hold on
        
        loc = x(locs(k));
        
        N = (find(x < (loc + dx/2) & x > (loc - dx/2)) - find(x==0));
        
        T_p = T_s(:,i_l+N : i_r+N);
        
        E_l2(j,k) = norm(T_i - T_p);
        %k
    end
    
    figure(4)
    contourf(x(i_l:i_r)+loc,z,T_i,15)
    hold on
    contour(x(i_l+N:i_r+N),z,T_p,15,'LineWidth', 1, 'Color', 'k')
    
    axis equal tight
    
    
    
end

figure(2)
for k = 1:size(locs,2)
    plot(t,E_l2(:,k),'Linewidth',2);
    hold on
end