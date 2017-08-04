clear all;
close all;

fileID = fopen('CloseFormulaGaussian.txt','r');
closedPrice = fscanf(fileID, '%f');
fclose(fileID);

MC_iterations = [100, 500, 1000, 5000, 10000, 50000, 100000];


col = ['k','r','b','m','g','c','y'];

MC_DOWN_0 = [];
MC_UP_0 = [];
MC_EST_0 = [];

for i=1:size(MC_iterations,2)
  file_name = strcat('MonteCarloGaussian',num2str(MC_iterations(i)),'Iterations.txt');
  fileID = fopen(file_name, 'r');
  res = fscanf(fileID, '%f');
  fclose(fileID);
  
   
  MC_DOWN_0 = [MC_DOWN_0, res(2) - res(4)];
  MC_UP_0 = [MC_UP_0, res(2) + res(4)];
  MC_EST_0 = [MC_EST_0, res(2)];
  
end


figure(1);
hold on;
p1 = plot([1:1:size(MC_iterations,2)], MC_DOWN_0, 'c');
p2 = plot([1:1:size(MC_iterations,2)], MC_UP_0, 'c');
p3 = plot([1:1:size(MC_iterations,2)], MC_EST_0, 'b');
p4 = plot([1:1:size(MC_iterations,2)], closedPrice*ones(1,size(MC_iterations,2)),'k');
hold off;
%title('theta = 0, M of MonteCarlo varies, BS Call');
legend([p3,p1,p4],'Estimated EL','IC' ,'Closed Formula', 'Location', 'NorthEast');
xlabel('M');

