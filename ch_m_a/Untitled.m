% Get x coords 
x1 = dlmread('D:\proga\ch_m_a\txe3.txt'); 

% Get y coords 
y1 = dlmread('D:\proga\ch_m_a\tye3.txt'); 
y2 = dlmread('D:\proga\ch_m_a\tye4.txt'); 
y3 = dlmread('D:\proga\ch_m_a\tye5.txt'); 

% cos 
xlabel('X'); 
ylabel('Y'); 
x = 0 : pi/100 :  pi * 3/10; 
y = tan(x); 
plot(x, y, '-k', 'LineWidth', 2) 
hold on; 

% Draw graphics 
% 5 
plot(x1, y1, 'b'), grid; 
hold on; 

plot(x1, y2, 'r'), grid; 
hold on; 

plot(x1, y3, 'g'), grid; 
hold on; 

legend('tg(x)', '3', '4', '5', 'Location','northwest');
