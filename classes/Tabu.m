classdef Tabu
    %TABU is class representing the tabu moves
    
    properties
        station %station cannot be moved
        route   %route where the station cannot be moved
        theta   %total duration of the tabu
        t       %number of iteration since the creation of the tabu 
    end
    
    
    methods
        function obj = Tabu(station,route,theta_min, theta_max)
            %TABU Construct an instance of this class
            if nargin==0
                return
            end
            obj.station =  station;
            obj.route = route;
            %duration of the tabu if uniform random integer between theta_min and theta_max 
            obj.theta = randi([theta_min,theta_max]);
            obj.t=0;
        end
        
        function obj = next(obj)
            % NEXT updates the tabu array after an iteration 
            for tabuMove=1:length(obj)
                obj(tabuMove).t = obj(tabuMove).t + 1;
            end
            still_tabu=[obj.t]<=[obj.theta];
            obj=obj(still_tabu);
        end

        function bool= is_tabu(obj,station,route)
            % tells if a move is tabu
            for tabu_move=obj
                if isequal(station,tabu_move.station) && route==tabu_move.route
                    bool = true;
                    return
                end
            end
            bool = false;
        end

    end
end

