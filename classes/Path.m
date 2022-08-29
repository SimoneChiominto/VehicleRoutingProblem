classdef Path < handle
    % path is a sequence of links in the complete graph of the stations in
    % the city
    properties
        stations        % array of class station
        n_stations      % number of stations in the path
        city            % objet of class city
    end

    methods

        function obj=Path(city,inds)
            % Path: constructor of objects of class Path
            % input:
            %   city: object of class city 
            %   inds: integer array, it represents the indices in the
            %       array of stations in the object city that are in the
            %       path
            if nargin==0
                return
            end
            if length(inds)>1
                % controllo se due indici consecutivi consecutivi sono
                % uguali. Se sono uguali allora cancello il secondo
                check= logical([1,inds(1:end-1)~=inds(2:end)]); 
                inds=inds(check);
            end
            obj.stations = city.stations(inds);
            obj.n_stations=length(inds);
            obj.city=city;
        end


        function paths=invert(obj)
            % invert: generate a paths in opposite direction

            % manage empty path array
            if isempty(obj)
                paths=obj;
                return
            end

            % manage path array
            paths(length(obj))=Path();  % preallocate path array array
            for i=1:length(obj)
                inds=[obj(i).stations.n];
                inds=inds(end:-1:1);
                paths(i)=Path(obj(i).city,inds);
            end
        end



        function path=new_concatenate(obj)
            % new_concatenate:  generate a new path concatenating paths in 
            %                   a path array
            
            % delete null paths
            null_path=[obj.n_stations]<=1;
            obj=obj(~null_path);

            starting_stations=zeros(1,length(obj)); % array of indices of starting stations
            ending_stations=zeros(1,length(obj));   % array of indices of ending stations
            for i=1:length(obj)
                starting_stations(i)=obj(i).stations(1).n;
                ending_stations(i)=obj(i).stations(end).n;
            end

            visited_inds=zeros(1,length(obj));      % boolean array of the paths already concatenated
            inverted_paths=zeros(1,length(obj));    % boolean array of the paths that are inverted to be concatenated
            inds=[];                                % indices of the stations in the concatenated path

            finish=min(starting_stations);
            for i=1:length(obj)
                start=finish;
                if ~any( starting_stations(~visited_inds)==start)
                    if ~any(ending_stations(~visited_inds)==start )
                        % se non trovo nè tra gli inizi nè tra le fini dei
                        % path una stazione con cui collegarsi ho un errore
                        error('non è possibile concatenare questi cammini')
                    else
                        ind=find(ending_stations==start & ~visited_inds,1);
                        inverted_paths(ind)=1;
                    end
                else
                    ind=find(starting_stations==start & ~visited_inds,1);
                end
                visited_inds(ind)=1;

                if inverted_paths(ind)==1
                    ns=[obj(ind).stations.n];
                    ns=ns(end:-1:1);
                    inds=[inds,ns];
                else
                    inds=[inds,obj(ind).stations.n];
                end
                finish=inds(end);
            end

            path=Path(obj(1).city,inds);
        end
        


        function paths=break_path(obj,breaking_points) 
            % break_paths:  break path in paths deleting the link after the
            %               point indicated
            %       input:
            %       obj:                object of class path
            %       breaking_points:    indices where to break the path
            last_break=0;
            breaking_points=sort(breaking_points);  % i punti di rottura devono essere ordinati
            paths(length(breaking_points)+1)=Path();    % preallocate path array
            for i=1:length(breaking_points)
                paths(i)=Path(obj.city,[obj.stations(last_break+1:breaking_points(i)).n]);
                last_break=breaking_points(i);
            end
            paths(length(breaking_points)+1)=Path(obj.city,[obj.stations(last_break+1:obj.n_stations).n]);
        end



        function routes=from_path_to_routes(obj)
            % from_path_to_routes: rende un path chiuso un oggetto di classe             
            %                       route    
           
            % manage empty array
            if isempty(obj)
                routes=route.empty;
                return
            end

            % manage route array
            %routes(length(obj))=route(); % preallocate output
            for i=1:length(obj)
               routes(i)=route(obj(i).city,[obj(i).stations.n]) ;
            end
        end



        function fig=plot(obj)
            % Plot: Visualization method
            fig=plot(obj(1).city);
            fig;
            hold on;
            for i=1:length(obj)
                plot([obj(i).stations.x],[obj(i).stations.y]);
            end
            hold off;
        end



        function paths=concatenate(obj,ords,directions)
            % variante di concatenate usata in 3 opt
            n=1;
            for i=1:size(ords)
                for j=1:size(directions)
                    appo_paths=obj;
                    to_be_inverted=find(directions(j,:));
                    appo_paths(to_be_inverted)=appo_paths(to_be_inverted).invert();
                    inds=[];
                    for k=1:4
                        inds=[inds,[appo_paths(k).stations.n]];
                    end
                    paths(n)=Path(obj(1).city,inds);
                    n=n+1;
                end
            end
        end

        function stations=nbhd(obj,station,p)
            route_stations=[obj.stations.n];
            if route_stations(1)==route_stations(end)
                route_stations=route_stations(1:end-1);
            end
            dist=obj.city.distances(station.n,route_stations);
            [~,inds]=sort(dist);
            len=length(inds);
            if dist(inds(1))==0 %significa che la stazione di cui vedere il vicinato è nella route
                stations=[obj.stations(inds(2:min([len,p+1])))];
                return
            end
            stations=[obj.stations(inds(1:min([len,p])))];
        end
    end
end