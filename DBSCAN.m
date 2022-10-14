function [class] = DBSCAN(data,MinPt)
%% Density-based spatial clustering of applications with noise (DBSCAN) 
% to cluster the candidates of parameter sets
%  given a set of points in some space, 
%  it groups together points that are closely packed together (points with many nearby neighbors), 
%  marking as outliers points that lie alone in low-density regions (whose nearest neighbors are too far away).

% Input:    
%    - data: data needs to be clustered. Here it should be the candidates of identified parameters.
%	 - MinPts: Minimum points considered as a cluster.  
%	 - Eps: a point a is directly reachable from point b if a is within the distance Eps

% output:
%     - class: indicate which class each point belongs
if nargin == 1;
	MinPts = 5;	
    Eps = sqrt((49.5^2*2+4^2));
end  
coef = 0.05;
[m,n] = size(data);  
x = [(1:m)' data];  
[m,n] = size(x);
types = zeros(1,m);
class = ones(1,m)*(-1);     % initialize all points as outliers
while isempty(find(class==1));
    coef = coef + 0.01;  
    Eps = Eps*coef;  % by default,EPS = sqrt((49.5^2*2+4^2))*0.06;
    dealed = zeros(m,1); 
    dis = calDistance(x(:,2:n));  % calculate the distance between any pairs of points
    number = 1;
    for i = 1:m  
    %
        if dealed(i) == 0   % if point i has not been dealed
            xTemp = x(i,:);  
            D = dis(i,:);
            ind = find(D<=Eps);
			%
            if length(ind)>1 && length(ind)<MinPts+1  
                types(i) = 0;  
                class(i) = 0;  
            end   
			%
            if length(ind) == 1  
                types(i) = -1;  
                class(i) = -1;  
                dealed(i) = 1;  
            end  
            %
            if length(ind)>=MinPts+1  
                types(xTemp(1,1)) = 1;  
                class(ind) = number;  
                % 
                while ~isempty(ind)  
                    yTemp = x(ind(1),:);  
                    dealed(ind(1)) = 1;  
                    ind(1) = [];  
                    D = dis(yTemp(1,1),:);%
                    ind_1 = find(D<=Eps);  
                    if length(ind_1)>1%
                        class(ind_1) = number;  
                        if length(ind_1) >= MinPts+1  
                            types(yTemp(1,1)) = 1;  
                        else  
                            types(yTemp(1,1)) = 0;  
                        end  

                        for j=1:length(ind_1)  
                           if dealed(ind_1(j)) == 0  
                              dealed(ind_1(j)) = 1;  
                              ind=[ind ind_1(j)];     
                              class(ind_1(j))=number;  
                           end                      
                       end  
                    end  
                end  
                number = number + 1;  
            end  
        end  
    end  
end
  
%% 
ind_2 = find(class==0);  
class(ind_2) = -1;  
types(ind_2) = -1;  
  
end
%%   
function [ dis ] = calDistance( x )  
    [m,n] = size(x);  
    dis = zeros(m,m);       
    for i = 1:m  
        for j = i:m  
            %  
            tmp =0;  
            for k = 1:n  
                tmp = tmp+(x(i,k)-x(j,k)).^2;  
            end  
            dis(i,j) = sqrt(tmp);  
            dis(j,i) = dis(i,j);  
        end  
    end  
end  
%% epsilon
function [Eps] = epsilon(x,k)  
  
% Function: [Eps] = epsilon(x,k)  
%  
% Aim:   
% Analytical way of estimating neighborhood radius for DBSCAN  
%  
% Input:   
% x - data matrix (m,n); m-objects, n-variables  
% k - number of objects in a neighborhood of an object  
% (minimal number of objects considered as a cluster)  
 
[m,n]=size(x);  
  
Eps=((prod(max(x)-min(x))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);  

end