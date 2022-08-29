classdef VRP
    %VRP is the class for the problem setting and solution of VRP. 
    
    properties
        city
        n_vehicles
        capacity
        method
        solution
    end
    
    methods
        function obj = VRP(city, n_vehicles, capacity) 
            % VRP class constructor  
            % input: 
            %   city : region of the stations of the vrp problem
            %   n_vehicles: number of maximum available vehicles
            %   capacity: maximum capacity of each vehicle

            if nargin==0
                % construct an empty vrp problem
                return
            end
            obj.n_vehicles = n_vehicles;
            obj.city = city;
            obj.capacity = capacity;  
        end


        function [obj,time]=solve(obj,method)
            % SOLVE: solving method for vrp problem
            % input:
            %   obj: vrp object to be solved
            %   method: string representing the name of the method the user
            %           want to use for solving the vrp or function handle
            %           of a custom function for solving vrp
            % output:
            %   obj: solved vrp
            %   time: total solving time

            init=tic;
            if nargin==1
                % default value for method is Clarke and Wright algorithm
                method='Clarke-Wright' 
            end

            if isa(method, 'function_handle')
                % user can give his/her solution method through a function
                % handle
                solving_function= method;
                obj.method = method;
            else

                switch method
                    case 'Clarke-Wright'
                        solving_function= @Clarke_Wright_Savings;
                        obj.method = @ Clarke_Wright_Savings;
                    case 'tabu-search'
                        solving_function= @TABUROUTE;
                        obj.method = @TABUROUTE;
                    otherwise
                        solving_function = @clusters;
                        obj.method = @clusters;
                end
            end
            obj=solving_function(obj);
            time=toc(init);
        end
     
        
        %
        % Solving methods
        %
        
        function obj = Clarke_Wright_Savings(obj)
            n_routes = obj.city.n_stations-1; % the real number of initial routes
            routes = route.empty;
            
            for i = 1:obj.city.n_stations 
                routes(end+1) = route(obj.city,[1,i,1]); % create routes (1,i,1)
                % the first route in routes is "[1,1,1]", we will delete it
                % at the end
            end
            
            % create savings matrix
            savings_matrix = zeros(obj.city.n_stations,obj.city.n_stations);
                for i = 1:obj.city.n_stations
                    for j = i+1:obj.city.n_stations
                        savings_matrix(i,j)=obj.city.distances(i,1)+obj.city.distances(1,j)-...
                            obj.city.distances(i,j);
                    end
                end
            savings_matrix(1,:)=zeros(1,obj.city.n_stations); % the depot is already in all the routes

            while n_routes > obj.n_vehicles 
                max_saving = max(savings_matrix(:)); 

                % if the while constraint is never met

                if max(savings_matrix(:))==0
                    %if no solution is found it tries with more vehicles
                    try
                        error('No solution found: number of vehicles too low');
                    catch
                        warning('No solution found: number of vehicles too low')
                        obj.n_vehicles=obj.n_vehicles+1;
                        obj=Clarke_Wright_Savings(obj);
                        return
                    end
                end


                
                % find largest saving
                [k,l] = find(savings_matrix == max_saving); 
                
                % in case there are two maxima
                if length(k)>1 && length(l)>1
                    k = k(1);
                    l = l(1);
                end
                
                %it will be useful when merging routes
                switched = 0;
                if routes(l).n_stations == 3
                    switched = 1;
                    temp = l;
                    l = k;
                    k = temp;
                end
                            
                
                % check if stations k and l are interior to a route
                % if neither of them is we can create a new provisional
                % route, else ignore this pair
                
                is_interior = 0;
                if (routes(k).stations(2).n ~= k) &&  (routes(k).stations(routes(k).n_stations-1).n ~= k)                    
                    is_interior = 1;
                elseif (routes(l).stations(2).n ~= l) &&  (routes(l).stations(routes(l).n_stations-1).n ~= l)
                    is_interior = 1;
                end
                
                % check if stations k and l are already into the same route
                % if not we can create a new provisional
                % route, else ignore this pair
                is_existent = 0;
                if routes(k).n_stations == routes(l).n_stations
                    routes_k_list = [routes(k).stations.n];
                    routes_l_list = [routes(l).stations.n];
                    if routes_k_list == routes_l_list
                        is_existent = 1;
                    end
                end
                
                
                if is_interior == 1 || is_existent ==1                  
                    if switched ==1
                        savings_matrix(l,k) = 0;
                    else
                        savings_matrix(k,l) = 0;
                    end
                    %go to the next pair
                    
                else
                    
                    merging_routes = [routes(k),routes(l)];
                    
                    % insert the new arc in the right place
                    if routes(k).stations(routes(k).n_stations-1).n==k && routes(l).stations(2).n==l
                        new_route = merging_routes.merge; 
                    elseif routes(k).stations(2).n==k && routes(l).stations(2).n==l
                        new_route = merging_routes.merge(2); 
                    elseif routes(k).stations(routes(k).n_stations-1).n==k && routes(l).stations(routes(l).n_stations-1).n==l
                        new_route = merging_routes.merge(3); 
                    else 
                        new_route = merging_routes.merge(4); 
                    end
                   
                                  
  
                    % check if the capacity constraint is met
                    if new_route.demand < obj.capacity
                    
                        index_array=[];
                        for i=1:new_route.n_stations-1 
                            index_array(end+1)=new_route.stations(i).n;
                        end
                
                        % associate each array element to the route in
                        % which its index is
                        for i = 1:obj.city.n_stations
                            if ismember(i,index_array)
                                routes(i)=new_route;
                            end
 
                        end
                    
                        n_routes = n_routes-1;
                    end
                    
                    % eliminate pair (k,l)
                    if switched ==1
                        savings_matrix(l,k) = 0;
                    else
                        savings_matrix(k,l) = 0;
                    end
                end
                
            end

            % eliminate double routes in the array
            old_routes = [];
            for i = 2:obj.city.n_stations
                for j = i+1:obj.city.n_stations
                    if routes(i)==routes(j)
                        if ~ismember(j,old_routes)
                            old_routes(end+1)=j;
                        end 
                    end
                end
            end
            routes(old_routes)=[];
            routes(1)=[];
            obj.solution = VrpSolution(routes);
        end



        function seeds = generate_seeds(obj) 
            %this function generates seed for clustering algorithm
            seeds_ind = [];
            max_dist = max(obj.city.distances(:,1));
            ind_max_dist = find(obj.city.distances(:,1)==max_dist);
            seeds_ind(end+1) = ind_max_dist;

            for k = 2:obj.n_vehicles
                dist = zeros(1,obj.city.n_stations);
                for j = 2:obj.city.n_stations
                    if ~ismember(j,seeds_ind)
                        dist(j) = obj.city.distances(1,j);
                        for i = 1:length(seeds_ind)
                            dist(j)=dist(j)+obj.city.distances(j,seeds_ind(i));
                        end
                    end
                end
                ind_max_dist = max(dist);
                seeds_ind(end+1)= find(dist==max(dist));
            end
            
            seeds = seeds_ind;
            
        end
        

        function [obj, seeds]=clusters(obj, seeds) 
            % CLUSTERS: this solving method use a two phases approach for solving
            % vrp. First it generates the clusters and then solve the tsp
            % problems quasi-optimally through the GENI algorithm 
            
            if nargin<2
                seeds = obj.generate_seeds; 
            end
           
            g = obj.city.distances([obj.city.stations.n],1)+obj.city.distances([obj.city.stations.n],seeds)-obj.city.distances(1,seeds); 

            problem = optimproblem('ObjectiveSense', 'min');
            isVisitedInTour = optimvar('isVisitedInTour', obj.city.n_stations, obj.n_vehicles, ...
                    'Type', 'integer', 'LowerBound', 0, 'UpperBound', 1);

            problem.Objective = sum(g.*isVisitedInTour , 'all');

            intoConstraint = optimconstr(obj.city.n_stations);

            intoConstraint(1)= sum(isVisitedInTour(1,:))==obj.n_vehicles;
            intoConstraint(2:end) = sum(isVisitedInTour(2:end,:),2)==1;
            problem.Constraints.intoConstraint = intoConstraint;
            %show(intoConstraint)


            seedConstraint = optimconstr(obj.n_vehicles,obj.n_vehicles);
            seedConstraint = isVisitedInTour(seeds,:)==eye(obj.n_vehicles);
            problem.Constraints.seedConstraint = seedConstraint;

            loadConstraint = optimconstr(obj.n_vehicles);

            for i = 1:obj.n_vehicles
                loadConstraint(i) = sum([obj.city.stations.demand]*isVisitedInTour(:,i)) <= obj.capacity;
            end 
            problem.Constraints.loadConstraint = loadConstraint;
            try
                [sol,minDist,exitflag,output] = solve(problem);
            

            matrix = round(sol.isVisitedInTour);

            routes = route.empty;
            for i=1:obj.n_vehicles
                stations_in_route = [find(matrix(:,i)==1)]';
                routes(end+1)=route(obj.city,[stations_in_route,1]);               
            end

            catch
                warning('No feasible solution found')
                obj.n_vehicles=obj.n_vehicles+1;
                obj=clusters(obj);
                return
            end

            
            %obj.solution = VrpSolution(routes);
                disp('optimizing solution')
            for k=1:length(routes)
                fprintf("\n optimizing route %d of %d",k,length(routes));
                if routes(k).n_stations>3
                   stations = [routes(k).stations];
                   routes(k)=GENI(stations(1:end-1),5);
                %routes(k)=US(routes(k),4);
                end

            end
            obj.solution = VrpSolution(routes);
        end
    end
end
     