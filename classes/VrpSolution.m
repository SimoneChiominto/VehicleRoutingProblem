classdef VrpSolution
    %VrpSolution is a class containing all information on a solution of a
    %vehicle routing problem
    
    properties
        routes      % array of the routes of the solution
        cost        % total cost of the solution
        n_vehicles  % number of veichles used
    end
    

    %%%%%%%%%%%%%%%%%%%
    %   METHODS
    %%%%%%%%%%%%%%%%%%%

    methods
        function obj = VrpSolution(routes)
            % VrpSolution Construct an instance of this class
            if nargin==0
                return
            end
            obj.routes=routes;
            obj.cost=sum([routes.cost]);
            obj.n_vehicles=length(routes);
        end

        
        function [Rr,indRr]=find_route(obj,station)
            %find_route finds the route in a vrpsolution where a station is located
            for r=1:length(obj.routes)
                if any(station.n==[obj.routes(r).stations.n])
                    Rr=obj.routes(r);
                    indRr=r;
                    return
                end
            end
            Rr=route.empty;
            indRr=r;
        end

        function bool=is_feasible(obj,capacity)
            %is_feasible tells if a solution is feasible given a maximum
            %capacity of the vehicles Q
            bool=all([obj.routes.demand]<=capacity);
        end
        
        function fig = plot(obj)
            %PLOT Visualization method
            fig = plot(obj.routes);
        end

        function out=copy(obj)
            % copies the values of a a vrp solution in another VrpSolution
            % objet. It is necessary since route is an handle class. 
            out=obj;
            for i=1:length(obj.routes)
                out.routes(i)=copy(obj.routes(i));
            end        
        end

    end
end

