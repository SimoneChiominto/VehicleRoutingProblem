classdef City  < handle
    % City is a array of stations with a matrix of distrances precomputed 
    properties
        stations    % array of stations                  
        distances   % matrix of the distances between stations
        n_stations  % numbero of stations in the city 
    end

    methods
        function obj=City(stations)
            % City: constructor of objects of class route
            % input:
            %   stations: array of object of class station
            if nargin==0
                return
            end
            % constructor
            obj.stations=stations;
            dist=zeros(length(stations));
            for i=1:length(stations)
                obj.stations(i).n=i;
                for j=i:length(stations)
                    dist(i,j)=((stations(i).x-stations(j).x)^2+...
                        (stations(i).y-stations(j).y)^2)^0.5;
                end
            end
            obj.distances=dist+dist';
            obj.n_stations=length(stations);
            for i=1:obj.n_stations
                obj.stations(i).city=obj;
            end
        end

        %
        % INITIALIZATION METHODS
        %

        function obj=initialize(~,North,East,Capacity)
            % create city from the coordinates
            len=length(North);
            if nargin<4
                Capacity=zeros(len,1);
            end
            
            stations(len)=station(); % allocate memory      
            stations(1)=station(North(1),East(1),Capacity(1),'depot');
            for i=2:len
                stations(i)=station(North(i),East(i),Capacity(i));
            end
            obj=City(stations);
        end

        
        function obj = initialize_from_file(obj, file)
            %transforms the file into an array of char
            fid = fopen(file);
            line1 = fgetl(fid);
            res=line1;
            while ischar(line1)
                line1 = fgetl(fid);
                res = char(res,line1);
            end
            fclose(fid);

            s = size(res);
            for i= 1:s(1)
                if strfind(res(i,:),'NODE_COORD_SECTION')==1
                    start_coor = i+1;
                elseif strfind(res(i,:),'DEMAND_SECTION')==1
                    end_coor = i-1;
                    start_demand = i+1;
                elseif strfind(res(i,:),'DEPOT_SECTION')==1
                    end_demand = i-1;
                    start_depot = i+1;
                end
            end

            coordinates = zeros(end_coor-start_coor+1,3);
            space = start_coor;
            for i = start_coor:end_coor
                coordinates(i-space+1,:) = str2num(res(i,:));
            end

            demand = zeros(end_demand-start_demand+1,2);
            space = start_demand;
            for i = start_demand:end_demand
                demand(i-space+1,:) = str2num(res(i,:));
            end

            depot = str2num(res(start_depot,:));
            if depot(1) ~= 1
                temp = coordinates(1,:);
                coordinates(1,:) = coordinates(depot(1),:);
                coordinates(depot(1),:) = coordinates(1,:);
            end

            obj = obj.initialize(coordinates(:,2),coordinates(:,3),demand(:,2));
            
        end
        
        %
        % OTHER METHODS
        %

        function obj= add_station(obj,station)
            % add_station: add a new station to the city.
            obj.n_stations = obj.n_stations +1;
            obj.stations(obj.n_stations) = station;
            dist=zeros(obj.n_stations-1,1);
            for i=1:obj.n_stations-1
                dist(i)=((obj.stations(i).x-station.x)^2+...
                    (obj.stations(i).y-station.y)^2)^0.5;
            end
            obj.distances=[obj.distances, dist;
                dist',0];
        end

        function fig=plot(obj)
            % simple visualization method
            fig=figure;
            hold on;
            names=string(1:obj.n_stations);
            sz=15;
            scatter([obj.stations.x],[obj.stations.y],sz,"filled");
            cond=string({obj.stations.category})=='depot';
            scatter([obj.stations(cond).x],[obj.stations(cond).y],sz,'red','filled');
            text([obj.stations.x]- 0.25, [obj.stations.y] + 0.25, names, 'HorizontalAlignment', 'right');
            axis equal;
            xlabel('East');
            ylabel('North');
            hold off
            return
        end
    end
end
  