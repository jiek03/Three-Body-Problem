clc, clear,close all
hold on
m1 = 1.99e30; %Sun mass
m2 = 1.898e27; %jovian mass
m3 = 7.35e18; %asteroid mass
s_a = 60^2*24*365; %years to seconds conversion factor
G = 6.67e-11 * s_a^2; %Universal Gravitation Constant
t1 = 0;
tf = 25; %años 
dt = 0.0001;
P = 1000; %Modulation of animation velocity
UA = 150e9;
theta = 3/4*pi - 1.2382; %Angular initial conditions of Jupiter
phi = 6.1; %Angular initial conditions of asteroid
%Function that calculates distance between two points
r_M = @(x1,y1,z1,x2,y2,z2) sqrt((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2);
%Acceleration Functions  
a_1 = {@(x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3) (-G*m2*(x_1-x_2)/r_M(x_1,y_1,z_1,x_2,y_2,z_2)^3 -G*m3*(x_1-x_3)/r_M(x_1,y_1,z_1,x_3,y_3,z_3)^3);
    @(x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3)  (-G*m2*(y_1-y_2)/r_M(x_1,y_1,z_1,x_2,y_2,z_2)^3 -G*m3*(y_1-y_3)/r_M(x_1,y_1,z_1,x_3,y_3,z_3)^3);
    @(x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3)(-G*m2*(z_1-z_2)/r_M(x_1,y_1,z_1,x_2,y_2,z_2)^3  -G*m3*(z_1-z_3)/r_M(x_1,y_1,z_1,x_3,y_3,z_3)^3)} ;

a_2 = {@(x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3) (-G*m1*(x_2-x_1)/r_M(x_1,y_1,z_1,x_2,y_2,z_2)^3 -G*m3*(x_2-x_3)/r_M(x_2,y_2,z_2,x_3,y_3,z_3)^3) ;
    @(x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3) (-G*m1*(y_2-y_1)/r_M(x_1,y_1,z_1,x_2,y_2,z_2)^3 -G*m3*(y_2-y_3)/r_M(x_2,y_2,z_2,x_3,y_3,z_3)^3);
    @(x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3) (-G*m1*(z_2-z_1)/r_M(x_1,y_1,z_1,x_2,y_2,z_2)^3 -G*m3*(z_2-z_3)/r_M(x_2,y_2,z_2,x_3,y_3,z_3)^3)};

a_3 = {@(x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3) (-G*m1*(x_3-x_1)/r_M(x_1,y_1,z_1,x_3,y_3,z_3)^3 - G*m2*(x_3-x_2)/r_M(x_2,y_2,z_2,x_3,y_3,z_3)^3);
    @(x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3) (-G*m1*(y_3-y_1)/r_M(x_1,y_1,z_1,x_3,y_3,z_3)^3 - G*m2*(y_3-y_2)/r_M(x_2,y_2,z_2,x_3,y_3,z_3)^3);
    @(x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3) (-G*m1*(z_3-z_1)/r_M(x_1,y_1,z_1,x_3,y_3,z_3)^3 - G*m2*(z_3-z_2)/r_M(x_2,y_2,z_2,x_3,y_3,z_3)^3)};

%Initial conditions
r_1i = [0,0,0];
r_2i = UA*5.4*[cos(theta),sin(theta),0];
r_3i = UA*7*[cos(phi),sin(phi),0];

v_1i = [0,0,0];
v_2i = s_a*13000*[-sin(theta),cos(theta),0];
v_3i = 2*UA*[-0.3,0.77,0];

%Velocity Verlet Results

[t,r1,r2,r3,v1,v2,v3] = V(a_1,a_2,a_3,0,tf,r_1i,r_2i,r_3i,v_1i,v_2i,v_3i,dt);


t_lo = length(t);

T = zeros(1,length(t));
U = zeros(1,length(t));
L = zeros(3,length(t));

p = m1*v1 + m2*v2 + m3*v3;

for i = 1:length(t)
L(:,i) = cross(r1(:,i),m1*v1(:,i)) + cross(r2(:,i),m2*v2(:,i)) + cross(r3(:,i),m3*v3(:,i)) ;

T(i) =  0.5*m1*dot(v1(:,i),v1(:,i)) + 0.5*m2*dot(v2(:,i),v2(:,i)) + 0.5*m3*dot(v3(:,i),v3(:,i));
U(i) =  -G*m1*m2/(r_M(r1(1,i),r1(2,i),r1(3,i),r2(1,i),r2(2,i),r2(3,i))) - G*m2*m3/(r_M(r2(1,i),r2(2,i),r2(3,i),r3(1,i),r3(2,i),r3(3,i))) - G*m1*m3/(r_M(r1(1,i),r1(2,i),r1(3,i),r3(1,i),r3(2,i),r3(3,i)))  ;

end
    
%Graphs and animations
for i = round((t1*(length(t)- 1)/tf))+1:P: length(t)
clf

hold on
grid on


tiledlayout(1,2)
set(gcf, 'Position', get(0, 'Screensize'));
nexttile
hold on
p2 = plot3(r2(1,i),r2(2,i),r2(3,i),'-o','Color' ,'black','MarkerSize',11 );
p2.MarkerFaceColor = '#D95319';

p3 = plot3(r3(1,i),r3(2,i),r3(3,i),'-o','Color' , 'black');
p3.MarkerFaceColor = '#707B7C';


xlim([r2(1,i)-UA,r2(1,i)+UA])
ylim([r2(2,i)-UA,r2(2,i)+UA])
xlabel('x [meters]')
ylabel('y [meters]')
zlabel('z [meters]')
title('Three Bodies on Gravitational Interaction')
grid on

nexttile
hold on
grid on
plot3(r1(1,:),r1(2,:),r1(3,:),'Color','#EDB120') %Sun
plot3(r2(1,:),r2(2,:),r2(3,:),'Color','#D35400') %Jupyter
plot3(r3(1,:),r3(2,:),r3(3,:),'Color','black') %Asteroid

p1 = plot3(r1(1,i),r1(2,i),r1(3,i),'-o','Color' , '#F39C12','MarkerSize',12);
p1.MarkerFaceColor = '#F1C40F';

p2 = plot3(r2(1,i),r2(2,i),r2(3,i),'-o','Color' ,'black','MarkerSize',8 );
p2.MarkerFaceColor = '#D95319';

p3 = plot3(r3(1,i),r3(2,i),r3(3,i),'-o','Color' , 'black');
p3.MarkerFaceColor = '#707B7C';

xlabel('x [meters]')
ylabel('y [meters]')
zlabel('z [meters]')
title('Sun - Jupiter - Asteroid system')
if min([(r1(1,:)), (r2(1,:))]) ~= max([(r1(1,:)), (r2(1,:))])
    xlim([min([(r1(1,:)), (r2(1,:)),(r3(1,:))]),max([(r1(1,:)), (r2(1,:)),(r3(1,:))])])
    
    
end

if min([(r1(2,:)), (r2(2,:))]) ~= max([(r1(2,:)), (r2(2,:))])
    ylim([min([(r1(2,:)), (r2(2,:)),(r3(2,:))]),max([(r1(2,:)), (r2(2,:)),(r3(2,:))])])

end

if min([(r1(3,:)), (r2(3,:))]) ~= max([(r1(3,:)), (r2(3,:))])
    zlim([min([(r1(3,:)), (r2(3,:)),(r3(3,:))]),max([(r1(3,:)), (r2(3,:)),(r3(3,:))])])
   
end

view(245,45)


pause(0.01)


end
figure
hold on 
plot(t,T)
plot(t,U)
plot(t,U+T)
title('System Energy')
xlabel('Time [time units]')
ylabel('Energy [energy units]')
legend ('Kinetic Energy', 'Potential Energy', 'Mechanical Energy')

%VELOCITY VERLET FUNCTION
function [t,r1,r2,r3,v1,v2,v3] = V(a_1,a_2,a_3,t_1,t_f,r_1i,r_2i,r_3i,v_1i,v_2i,v_3i,h)
t = t_1:h:t_f;
r1 = zeros(3,length(t));
r2 = zeros(3,length(t));
r3 = zeros(3,length(t));
v1 = zeros(3,length(t));
v2 = zeros(3,length(t));
v3 = zeros(3,length(t));
v_m1 = zeros(3,length(t));
v_m2 = zeros(3,length(t));
v_m3 = zeros(3,length(t));

v_m1(1,1) = v_1i(1) + 0.5*h*a_1{1}(r_1i(1),r_1i(2),r_1i(3),r_2i(1),r_2i(2),r_2i(3),r_3i(1),r_3i(2),r_3i(3));
v_m1(2,1) = v_1i(2) + 0.5*h*a_1{2}(r_1i(1),r_1i(2),r_1i(3),r_2i(1),r_2i(2),r_2i(3),r_3i(1),r_3i(2),r_3i(3));
v_m1(3,1) = v_1i(3) + 0.5*h*a_1{3}(r_1i(1),r_1i(2),r_1i(3),r_2i(1),r_2i(2),r_2i(3),r_3i(1),r_3i(2),r_3i(3));

v_m2(1,1) = v_2i(1) + 0.5*h*a_2{1}(r_1i(1),r_1i(2),r_1i(3),r_2i(1),r_2i(2),r_2i(3),r_3i(1),r_3i(2),r_3i(3));
v_m2(2,1) = v_2i(2) + 0.5*h*a_2{2}(r_1i(1),r_1i(2),r_1i(3),r_2i(1),r_2i(2),r_2i(3),r_3i(1),r_3i(2),r_3i(3));
v_m2(3,1) = v_2i(3) + 0.5*h*a_2{3}(r_1i(1),r_1i(2),r_1i(3),r_2i(1),r_2i(2),r_2i(3),r_3i(1),r_3i(2),r_3i(3));

v_m3(1,1) = v_3i(1) + 0.5*h*a_3{1}(r_1i(1),r_1i(2),r_1i(3),r_2i(1),r_2i(2),r_2i(3),r_3i(1),r_3i(2),r_3i(3));
v_m3(2,1) = v_3i(2) + 0.5*h*a_3{2}(r_1i(1),r_1i(2),r_1i(3),r_2i(1),r_2i(2),r_2i(3),r_3i(1),r_3i(2),r_3i(3));
v_m3(3,1) = v_3i(3) + 0.5*h*a_3{3}(r_1i(1),r_1i(2),r_1i(3),r_2i(1),r_2i(2),r_2i(3),r_3i(1),r_3i(2),r_3i(3));

r1(:,1) = r_1i;
r2(:,1) = r_2i;
r3(:,1) = r_3i;

v1(:,1) = v_1i;
v2(:,1) = v_2i;
v3(:,1) = v_3i;


for i = 1:length(t) - 1 
    if i == 111500  %event
        v3(1,i) = 0.75*(v3(1,i-1)*cos(-pi/6) - v3(2,i-1)*sin(-pi/6 ));
        v3(2,i) = 0.95*(v3(2,i-1)*cos(-pi/6) + v3(1,i-1)*sin(-pi/6 ));
    end
    v_m1(1,i+1) = v1(1,i) + 0.5*h*a_1{1}(r1(1,i),r1(2,i),r1(3,i),r2(1,i),r2(2,i),r2(3,i),r3(1,i),r3(2,i),r3(3,i));
    v_m1(2,i+1) = v1(2,i) + 0.5*h*a_1{2}(r1(1,i),r1(2,i),r1(3,i),r2(1,i),r2(2,i),r2(3,i),r3(1,i),r3(2,i),r3(3,i));
    v_m1(3,i+1) = v1(3,i) + 0.5*h*a_1{3}(r1(1,i),r1(2,i),r1(3,i),r2(1,i),r2(2,i),r2(3,i),r3(1,i),r3(2,i),r3(3,i));
    
    v_m2(1,i+1) = v2(1,i) + 0.5*h*a_2{1}(r1(1,i),r1(2,i),r1(3,i),r2(1,i),r2(2,i),r2(3,i),r3(1,i),r3(2,i),r3(3,i));
    v_m2(2,i+1) = v2(2,i) + 0.5*h*a_2{2}(r1(1,i),r1(2,i),r1(3,i),r2(1,i),r2(2,i),r2(3,i),r3(1,i),r3(2,i),r3(3,i));
    v_m2(3,i+1) = v2(3,i) + 0.5*h*a_2{3}(r1(1,i),r1(2,i),r1(3,i),r2(1,i),r2(2,i),r2(3,i),r3(1,i),r3(2,i),r3(3,i));
    
    v_m3(1,i+1) = v3(1,i) + 0.5*h*a_3{1}(r1(1,i),r1(2,i),r1(3,i),r2(1,i),r2(2,i),r2(3,i),r3(1,i),r3(2,i),r3(3,i));
    v_m3(2,i+1) = v3(2,i) + 0.5*h*a_3{2}(r1(1,i),r1(2,i),r1(3,i),r2(1,i),r2(2,i),r2(3,i),r3(1,i),r3(2,i),r3(3,i));
    v_m3(3,i+1) = v3(3,i) + 0.5*h*a_3{3}(r1(1,i),r1(2,i),r1(3,i),r2(1,i),r2(2,i),r2(3,i),r3(1,i),r3(2,i),r3(3,i));
    
    r1(1,i+1) = r1(1,i) + h*v_m1(1,i+1);
    r1(2,i+1) = r1(2,i) + h*v_m1(2,i+1);
    r1(3,i+1) = r1(3,i) + h*v_m1(3,i+1);
    
    r2(1,i+1) = r2(1,i) + h*v_m2(1,i+1);
    r2(2,i+1) = r2(2,i) + h*v_m2(2,i+1);
    r2(3,i+1) = r2(3,i) + h*v_m2(3,i+1);
    
    r3(1,i+1) = r3(1,i) + h*v_m3(1,i+1);
    r3(2,i+1) = r3(2,i) + h*v_m3(2,i+1);
    r3(3,i+1) = r3(3,i) + h*v_m3(3,i+1);
    
    v1(1,i+1) = v_m1(1,i+1) + 0.5*h*a_1{1}(r1(1,i+1),r1(2,i+1),r1(3,i+1),r2(1,i+1),r2(2,i+1),r2(3,i+1),r3(1,i+1),r3(2,i+1),r3(3,i+1));
    v1(2,i+1) = v_m1(2,i+1) + 0.5*h*a_1{2}(r1(1,i+1),r1(2,i+1),r1(3,i+1),r2(1,i+1),r2(2,i+1),r2(3,i+1),r3(1,i+1),r3(2,i+1),r3(3,i+1));
    v1(3,i+1) = v_m1(3,i+1) + 0.5*h*a_1{3}(r1(1,i+1),r1(2,i+1),r1(3,i+1),r2(1,i+1),r2(2,i+1),r2(3,i+1),r3(1,i+1),r3(2,i+1),r3(3,i+1));
    
    v2(1,i+1) = v_m2(1,i+1) + 0.5*h*a_2{1}(r1(1,i+1),r1(2,i+1),r1(3,i+1),r2(1,i+1),r2(2,i+1),r2(3,i+1),r3(1,i+1),r3(2,i+1),r3(3,i+1));
    v2(2,i+1) = v_m2(2,i+1) + 0.5*h*a_2{2}(r1(1,i+1),r1(2,i+1),r1(3,i+1),r2(1,i+1),r2(2,i+1),r2(3,i+1),r3(1,i+1),r3(2,i+1),r3(3,i+1));
    v2(3,i+1) = v_m2(3,i+1) + 0.5*h*a_2{3}(r1(1,i+1),r1(2,i+1),r1(3,i+1),r2(1,i+1),r2(2,i+1),r2(3,i+1),r3(1,i+1),r3(2,i+1),r3(3,i+1));
    
    v3(1,i+1) = v_m3(1,i+1) + 0.5*h*a_3{1}(r1(1,i+1),r1(2,i+1),r1(3,i+1),r2(1,i+1),r2(2,i+1),r2(3,i+1),r3(1,i+1),r3(2,i+1),r3(3,i+1));
    v3(2,i+1) = v_m3(2,i+1) + 0.5*h*a_3{2}(r1(1,i+1),r1(2,i+1),r1(3,i+1),r2(1,i+1),r2(2,i+1),r2(3,i+1),r3(1,i+1),r3(2,i+1),r3(3,i+1));
    v3(3,i+1) = v_m3(3,i+1) + 0.5*h*a_3{3}(r1(1,i+1),r1(2,i+1),r1(3,i+1),r2(1,i+1),r2(2,i+1),r2(3,i+1),r3(1,i+1),r3(2,i+1),r3(3,i+1));
end
end
