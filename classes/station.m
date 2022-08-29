classdef station
    %STATION is a class representing a station to be visited in a tsp or
    %vrp
    
    properties
        x   
        y
        demand
        category
        city
        n
    end
    
    methods
        function obj = station(x,y,demand,category)
            %STATION Construct an instance of this class
            if nargin==0
                return
            end
            obj.x=x;
            obj.y=y;
            obj.demand=demand;
            if nargin>3
                obj.category=category;
                return
            end
            obj.category='station';
            
        end
       

    end
end

