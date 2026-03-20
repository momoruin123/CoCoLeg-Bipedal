function traj_vector = trajToVector(config, traj)
% 把轨迹扩展成整个stride的,前一半不动,后一半交换左右腿
x1 = traj.x;
u1 = traj.u;
t1 = traj.t;

% Data Phase 2
end_t = t1(end);
t2 = end_t + t1;

x2 = x1;
x2(:,2) = x1(:,3);
x2(:,3) = x1(:,2);
x2(:,4) = x1(:,5);
x2(:,5) = x1(:,4);
x2(:,7) = x1(:,8);
x2(:,8) = x1(:,7);
x2(:,9) = x1(:,10);
x2(:,10) = x1(:,9);

u2 = u1;
u2(:,1) = u1(:,2);
u2(:,2) = u1(:,1);
u2(:,3) = u1(:,4);
u2(:,4) = u1(:,3);

traj_vector.x = [x1;x2];
traj_vector.u = [u1;u2];
traj_vector.t = [t1;t2];
end