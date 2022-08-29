function vrp = TABUROUTE(vrp)
%Heuristic for solving VRP based on tabu search
city=vrp.city;

alpha=1;
F1_min=inf;
F2_min=inf;
theta_min=5;
theta_max=10;
h=10;
Q=vrp.capacity;

F1= @(sol) sol.cost;
F2= @(sol,alpha) F1(sol) +...
    alpha * sum ( max([sol.routes.demand]-Q,0) ); 


%
% Initialization
%
disp('Started computing initial solution')
tsp_sol=GENI(city.stations,5);

%tsp_sol=US(tsp_sol,5);

vrp_sol=naive_clusters(tsp_sol,vrp.capacity,vrp.n_vehicles);

disp('Finished computing initial solution')
if vrp_sol.is_feasible(vrp.capacity) && F1(vrp_sol)<F1_min
    F1_min=F1(vrp_sol);
    S_min=tsp_sol;
end

F2_min=min(F2(vrp_sol,alpha),F2_min);
if F2_min==F2(vrp_sol,alpha)
    S_tilde_min=tsp_sol;
end
vrp.solution=vrp_sol;





%
% First round of SEARCH algorithm
%

W=city.stations(2:end);
fv=zeros(1,length(city.stations));
q=2*vrp.n_vehicles;
p1=3;
p2=p1;
n_max=city.n_stations;
g=0.01;


disp('Started first use of SEARCH algorithm')
[sol_vrp,sol_tilde_vrp,F1_min,alpha,fv]=SEARCH(W,q,p1,p2,theta_min,theta_max,g,h,n_max...
   ,vrp,alpha,fv);
if isfinite(F1_min)
   S=sol_vrp;
else
   S=sol_tilde_vrp;
end
vrp.solution=S;

save('partial_test') ;


%
% Second round of SEARCH algorithm
%
n_max=5*city.n_stations;

disp('Started second use of SEARCH algorithm')
[sol_vrp,sol_tilde_vrp,F1_min,alpha,fv]=SEARCH(W,q,p1,p2,theta_min,theta_max,g,h,n_max...
   ,vrp,alpha,fv);
if isfinite(F1_min)
   S=sol_vrp;
else
   S=sol_tilde_vrp;
end
vrp.solution=S;


%
% Third round of SEARCH algorithm
%

[~,inds]=sort(fv);
inds=inds(end:-1:floor(city.n_stations/2));
W=city.stations(inds);
q=floor(city.n_stations/4);
n_max=city.n_stations;

disp('Started third use of SEARCH algorithm')

[sol_vrp,~,F1_min,~]=SEARCH(W,q,p1,p2,theta_min,theta_max,g,h,n_max...
    ,vrp,alpha,fv);


if isfinite(F1_min)
    disp('Started optimizing routes using US algorithm');
    for r=1:length(sol_vrp.routes)  
        fprintf('\nOptimizing route %d',r);
        sol_vrp.routes(r)=US(sol_vrp.routes(r),5); 
    end
    vrp.solution=sol_vrp;
    return
else
    error('no solution found')
end

end



function vrp_sol=naive_clusters(tsp_sol,capacity,m_max)
%
% Cluster Heuristic based of geometrical intuition.
%
clusters={};

total_demand=0;
n_clusters=1;
starting_index=2;

for i=2:tsp_sol.n_stations-1
    %se non posso aggiungere altri veicoli stop
    if n_clusters==m_max || i==tsp_sol.n_stations-1
        clusters{n_clusters}=tsp_sol.stations(starting_index:end);
        break
    end

    total_demand=total_demand+tsp_sol.stations(i).demand;
    
    if total_demand>capacity
        total_demand=tsp_sol.stations(i).demand;
        clusters{n_clusters}=tsp_sol.stations(starting_index:i-1);
        starting_index=i;
        n_clusters=n_clusters+1;
        
    end

end    

n_clusters=length(clusters);
routes(n_clusters)=route();
for i=1:n_clusters
    routes(i)=route( tsp_sol.city , [1,[clusters{i}.n],1] );
end

vrp_sol=VrpSolution(routes);

end