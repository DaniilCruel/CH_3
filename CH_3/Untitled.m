filename = 'Eiler_Implicit.txt';
[B, delimiterOut]= importdata(filename);
hold on;
grid on;
for i =2:4
    plot (B(:,1),B(:, i));
end