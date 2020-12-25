%% Parameters

number_states=13;
number_measurements=10;

alpha=[pi/4,0.523599,0.561489,0.659058,pi/4,0.261825,0.362367,0.561489,pi/4,0,0.261825,0.523599,pi/4];
theta=[pi/2,pi/2,0.24+pi/2,0.43+pi/2,0.47+pi/2,pi/2,0.65+pi/2,0.989+pi/2,0.9817+pi/2,0,0,0,0];

%% Analytical values

alpha_matrix=2.*alpha;
theta_matrix=2.*theta;

unity=[1 0 ; 0 1];
rho_i_th=zeros(4,4,number_states);
rho_A_i_th=zeros(2,2,number_states);

for i=1:13
        Ry=[cos(alpha(i)) -sin(alpha(i)) ; sin(alpha(i)) cos(alpha(i))];
        CRy=[1 0 0 0 ; 0 1 0 0 ; 0 0 cos(theta(i)) -sin(theta(i)); 0 0 sin(theta(i)) cos(theta(i))];
        U_prep=CRy*kron(Ry,unity);
        rho_i_th(:,:,i)=(U_prep*[1;0;0;0])*(U_prep*[1;0;0;0])';
        rho_A_i_th(:,:,i)=PartialTrace(rho_i_th(:,:,i),[2],[2,2]);
end

C=[0,0,0,0,0,0,0,0,0,0,0,0,0]';
VA=[0,0,0,0,0,0,0,0,0,0,0,0,0]';
PA=[0,0,0,0,0,0,0,0,0,0,0,0,0]';
for i=1:13
C(i)=concurrence(rho_i_th(:,:,i));
VA(i)=2*abs(rho_A_i_th(1,2,i));
PA(i)=abs(rho_A_i_th(1,1,i)-rho_A_i_th(2,2,i));
end

%% Read experimental results

% Initialization
tab_rho=zeros(4,4,number_states,number_measurements);
tab_rho_A=zeros(2,2,number_states,number_measurements);
tab_rho_B=zeros(2,2,number_states,number_measurements);
tab_lambda=zeros(number_states,number_measurements,4);
tab_VA=zeros(number_states,number_measurements);
tab_PA=zeros(number_states,number_measurements);
tab_C=zeros(number_states,number_measurements);
tab_JB=zeros(number_states,number_measurements);
mean_VA=zeros(number_states,1);
mean_PA=zeros(number_states,1);
mean_C=zeros(number_states,1);
mean_JB=zeros(number_states,1);
tab_rho_mitigated=zeros(4,4,number_states,number_measurements);
tab_rho_A_mitigated=zeros(2,2,number_states,number_measurements);
tab_rho_B_mitigated=zeros(2,2,number_states,number_measurements);
tab_lambda_mitigated=zeros(number_states,number_measurements,4);
tab_VA_mitigated=zeros(number_states,number_measurements);
tab_PA_mitigated=zeros(number_states,number_measurements);
tab_C_mitigated=zeros(number_states,number_measurements);
tab_JB_mitigated=zeros(number_states,number_measurements);
mean_VA_mitigated=zeros(number_states,1);
mean_PA_mitigated=zeros(number_states,1);
mean_C_mitigated=zeros(number_states,1);

% Raw values
for j=1:number_states
j=14-j; 
theta_i_str=num2str(14-j);

    for f=1:number_measurements
        f_str=num2str(f-1);
        measured_state = join(['Quantum_VPC_measurement_',f_str,'_state_',theta_i_str,'.txt']);
        vpc_file=fopen(measured_state,'r');
        vpc=str2num(fgets(vpc_file));
        vpc = load(measured_state);
        tab_rho(1,1,j,f)=vpc(:,8);
        tab_rho(1,2,j,f)=vpc(:,9);
        tab_rho(1,3,j,f)=vpc(:,10);
        tab_rho(1,4,j,f)=vpc(:,11);
        tab_rho(2,1,j,f)=vpc(:,12);
        tab_rho(2,2,j,f)=vpc(:,13);
        tab_rho(2,3,j,f)=vpc(:,14);
        tab_rho(2,4,j,f)=vpc(:,15);
        tab_rho(3,1,j,f)=vpc(:,16);
        tab_rho(3,2,j,f)=vpc(:,17);
        tab_rho(3,3,j,f)=vpc(:,18);
        tab_rho(3,4,j,f)=vpc(:,19);
        tab_rho(4,1,j,f)=vpc(:,20);
        tab_rho(4,2,j,f)=vpc(:,21);
        tab_rho(4,3,j,f)=vpc(:,22);
        tab_rho(4,4,j,f)=vpc(:,23);
        tab_rho_A(:,:,j,f)=PartialTrace(tab_rho(:,:,j,f),[2],[2,2]);
        tab_rho_B(:,:,j,f)=PartialTrace(tab_rho(:,:,j,f),[1],[2,2]);
        tab_VA(j,f)=2*abs(tab_rho_A(1,2,j,f));
        tab_PA(j,f)=abs((tab_rho_A(1,1,j,f)-tab_rho_A(2,2,j,f)));
        tab_C(j,f)=concurrence(tab_rho(:,:,j,f));
        tab_JB(j,f)=tab_C(j,f)*tab_C(j,f)+tab_VA(j,f)*tab_VA(j,f)+tab_PA(j,f)*tab_PA(j,f);
        tab_lambda(j,f,:)=sort(eig(tab_rho(:,:,j,f)));
        tab_JB(j,f)=(tab_C(j,f)*tab_C(j,f)+tab_VA(j,f)*tab_VA(j,f)+tab_PA(j,f)*tab_PA(j,f));
    
        fclose(vpc_file);

    end
    mean_VA(j)=mean(tab_VA(j,:));
    mean_PA(j)=mean(tab_PA(j,:));
    mean_C(j)=mean(tab_C(j,:));

end

% With qiskit ignis error mitigation
for j=1:number_states
j=14-j; 

theta_i_str=num2str(14-j);

    for f=1:number_measurements
        f_str=num2str(f-1);
        measured_state = join(['Quantum_VPC_measurement_mitigated_errors_',f_str,'_state_',theta_i_str,'.txt']);
        vpc = load(measured_state);
        tab_rho_mitigated(1,1,j,f)=vpc(:,8);
        tab_rho_mitigated(1,2,j,f)=vpc(:,9);
        tab_rho_mitigated(1,3,j,f)=vpc(:,10);
        tab_rho_mitigated(1,4,j,f)=vpc(:,11);
        tab_rho_mitigated(2,1,j,f)=vpc(:,12);
        tab_rho_mitigated(2,2,j,f)=vpc(:,13);
        tab_rho_mitigated(2,3,j,f)=vpc(:,14);
        tab_rho_mitigated(2,4,j,f)=vpc(:,15);
        tab_rho_mitigated(3,1,j,f)=vpc(:,16);
        tab_rho_mitigated(3,2,j,f)=vpc(:,17);
        tab_rho_mitigated(3,3,j,f)=vpc(:,18);
        tab_rho_mitigated(3,4,j,f)=vpc(:,19);
        tab_rho_mitigated(4,1,j,f)=vpc(:,20);
        tab_rho_mitigated(4,2,j,f)=vpc(:,21);
        tab_rho_mitigated(4,3,j,f)=vpc(:,22);
        tab_rho_mitigated(4,4,j,f)=vpc(:,23);
        tab_rho_A_mitigated(:,:,j,f)=PartialTrace(tab_rho_mitigated(:,:,j,f),[2],[2,2]);
        tab_rho_B_mitigated(:,:,j,f)=PartialTrace(tab_rho_mitigated(:,:,j,f),[1],[2,2]);
        tab_VA_mitigated(j,f)=2*abs(tab_rho_A_mitigated(1,2,j,f));
        tab_PA_mitigated(j,f)=abs((tab_rho_A_mitigated(1,1,j,f)-tab_rho_A_mitigated(2,2,j,f)));
        tab_C_mitigated(j,f)=concurrence(tab_rho_mitigated(:,:,j,f));
        tab_JB_mitigated(j,f)=tab_C_mitigated(j,f)*tab_C_mitigated(j,f)+tab_VA_mitigated(j,f)*tab_VA_mitigated(j,f)+tab_PA_mitigated(j,f)*tab_PA_mitigated(j,f);
        tab_lambda_mitigated(j,f,:)=sort(eig(tab_rho_mitigated(:,:,j,f)));

    end
    mean_VA_mitigated(j)=mean(tab_VA_mitigated(j,:));
    mean_PA_mitigated(j)=mean(tab_PA_mitigated(j,:));
    mean_C_mitigated(j)=mean(tab_C_mitigated(j,:));

end

%% Plots the Jakob-Bergou inequality in the V,P,C space

figure
for j=1:number_states
    scatter3(tab_VA(j,:),tab_PA(j,:),tab_C(j,:),'o','MarkerEdgeColor','black','MarkerFaceColor','white');
    hold on
end
axis equal
axis([0 1 0 1 0 1])
hold on
[x,y,z] = sphere(150);
hSurface=surf(x,y,z);
set(hSurface,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.3,'EdgeAlpha', 0.2);
xlabel('V_A')
ylabel('P_A')
zlabel('C')
title('Raw values')
view(135,20)

figure
for j=1:number_states
    scatter3(tab_VA_mitigated(j,:),tab_PA_mitigated(j,:),tab_C_mitigated(j,:),'o','MarkerEdgeColor','black','MarkerFaceColor','white');
    hold on
end
axis equal
axis([0 1 0 1 0 1])
hold on
[x,y,z] = sphere(150);
hSurface=surf(x,y,z);
set(hSurface,'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.3,'EdgeAlpha', 0.2);
xlabel('V_A')
ylabel('P_A')
zlabel('C')
title('With qiskit ignis error mitigation')
view(135,20)