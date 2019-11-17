%% Read data

% Number of points (number of different couples alpha, theta)
nb_points=13;
% Number of measurements (of V_K,P_k,C for each point)
nb_measurements=10;

% Initialize data vectors
alpha=zeros(nb_points*nb_measurements,1);
theta=zeros(nb_points*nb_measurements,1);
V_A=zeros(nb_points*nb_measurements,1);
P_A=zeros(nb_points*nb_measurements,1);
V_B=zeros(nb_points*nb_measurements,1);
P_B=zeros(nb_points*nb_measurements,1);
C=zeros(nb_points*nb_measurements,1);

% Read data
i=0;
for j=1:nb_measurements
    j_str=num2str(j-1);
for f=1:nb_points
    i=i+1;
    f_str=num2str(f);
    a_charger = join(['Quantum_VPC_measurement_',j_str,'_point_',f_str,'.txt']);
    vdc = load(a_charger);
    alpha(i)=vdc(:,1);
    theta(i)=vdc(:,2);
    V_A(i)=vdc(:,3);
    P_A(i)=vdc(:,4);
    V_B(i)=vdc(:,5);
    P_B(i)=vdc(:,6);
    C(i)=vdc(:,7);
end
end
%% Plot data

% Plot on V, P, C sphere for qubit A
figure
scatter3(V_A,P_A,C,'o','MarkerEdgeColor','black','MarkerFaceColor','red');
hold on
axis equal
axis([0 1 0 1 0 1])
[x,y,z] = sphere(150);
hSurface=surf(x,y,z);
set(hSurface,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.3,'EdgeAlpha', 0.2);
xlabel('V_A')
ylabel('P_A')
zlabel('C')
title('Qubit A')
view(135,20)