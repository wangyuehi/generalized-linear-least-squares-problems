
x = 20:20:500;


% hold on;
%  
% for i = 1:3
%     plot(x,time(:,i))
% end
% hold off

hold on;
 
for i = 1:3
    plot(x,relaError(:,i))
end
hold off