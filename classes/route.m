classdef route < Path 
    % ROUTE is an hamiltonyan cycle on a complete graph that starts and
    % ends with a depot

    properties
        demand  % demand is the total demand of the route
        cost    % cost is the total length of the route
    end


    %
    % METHODS
    %

    methods
        function obj = route(city,inds)
            % route: constructor of objects of class route
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

            obj.demand = sum([city.stations(inds).demand]);
            obj.cost =  sum(diag(city.distances([inds(1:end-1)],[inds(2:end)])));

            % check per controllare che l'inizio e la fine del ciclo siano
            % lo stesso e siano un depot
            
            check=(isequal(obj.stations(1),obj.stations(end)));
            if check==0
                error("Starting station is not a depot or it is not a cycle")
            end
        end


        function obj=invert(obj)
            % gira una route nel senso inverso
            obj=invert@Path(obj);
            obj=obj.from_path_to_routes();
        end


        function ind=plus(obj,station,i)
            % addition in the additive group defined by the route
            ind=mod(station-1+i,obj.n_stations-1)+1;
        end

        function route = merge(obj, type)
            % merge:  generates a new route merging routes in 
            %                   a route array
            % INPUTS:
            % type = integer that indicates how many routes have to be
            %        inverted before merging
            %        1 -> merge (default)
            %        2 -> invert first route, then merge
            %        3 -> invert second route, then merge
            %        4 -> invert both routes, then merge
            
            stat_to_delete = obj(1).n_stations;
            demand_to_delete = obj(1).stations(stat_to_delete).demand;
            ind_stat_to_delete = obj(1).stations(stat_to_delete).n;
            
            if nargin <2
                type = 1;
            end
            
            if type == 2
                obj(1) = obj(1).invert;
            elseif type == 3
                obj(2) = obj(2).invert;
            elseif type == 4
                obj(1) = obj(1).invert;
                obj(2) = obj(2).invert;
            end
           
            route_array = obj.new_concatenate;
            new_route = route_array.from_path_to_routes;
            new_route.stations(stat_to_delete)=[];
            new_route.n_stations = new_route.n_stations-1;
            new_route.cost = new_route.cost -  new_route.city.distances(new_route.stations(stat_to_delete-1).n,ind_stat_to_delete)-...
                new_route.city.distances(new_route.stations(stat_to_delete).n,ind_stat_to_delete)+...
                new_route.city.distances(new_route.stations(stat_to_delete-1).n,new_route.stations(stat_to_delete).n);
            new_route.demand = new_route.demand-demand_to_delete;
            route = new_route;
                
        end

        



        %%%%%%%%%%%%%%%%%%%
        % Insertion methods
        %%%%%%%%%%%%%%%%%%%

        function obj = insert(obj,station,position)
            % insert:   methods to insert a station in a position in the
            %           route
            if nargin==2
                obj=obj.greedy_insert(station);
            end
            stations=[obj.stations(1:position-1),station,obj.stations(position:end)];
            %city=station.city;
            
            obj.stations=stations;
            obj.n_stations = obj.n_stations+1;

            dist=obj.city.distances;
            obj.cost=obj.cost +...
                dist(station.n,obj.stations(position).n) + ...
                dist(station.n,obj.stations(position-1).n) - ...
                dist(obj.stations(position).n,obj.stations(position-1).n);
            
            obj.demand=obj.demand+station.demand;

            %obj=route(city,[stations.n]);
        end


        function obj = append(obj,station)
            % append: it is equaivalent to insert in the last position
            obj=obj.insert(station,obj.n_stations);
        end


        function obj=greedy_insert(obj,station)
            % greedy_insert: find the place where extramiliage is minimized
            %                and use insert
            dist=obj.city.distances;
            extra_miles= -diag(dist([obj.stations(1:end-1).n],[obj.stations(2:end).n]))...
                + dist([obj.stations(1:end-1).n],station.n)...
                + dist([obj.stations(2:end).n],station.n);
            ind=find(extra_miles==min(extra_miles),1)+1;
            obj=obj.insert(station,ind);
        end


        function obj=type1_insertion(obj,station,i,j,k)
            % type 1 insertion is described in the article by M.Gendreau,
            % A.Hertz,G.Laporte 'New insertion and post-optimization procedures
            % for the traveling salesman problem.' 
            
            %notare che gli indici i,j,k sono gli indici della route non
            %della città
            cond= k~=i && k~=j && i~=j && all(station.n~=[obj.stations.n]);
            if ~cond
                error("Type I insertion conditions not verified");
            end
            %
            len=obj.n_stations-1;
            if obj.is_in_the_middle(i,j,k)
                obj=obj.invert();
                inds=mod(len-[i-1,j-1,k-1],len)+1;
                i=inds(1);
                j=inds(2);
                k=inds(3);
            end

            iplus1=obj.plus(i,1);
            jplus1=obj.plus(j,1);
            kplus1=obj.plus(k,1);

            paths=obj.break_path([i,j,k]);
            new_links=Path(obj.city,[obj.stations(i).n,station.n]);
            new_links(2)=Path(obj.city,[station.n ,obj.stations(j).n]);
            new_links(3)=Path(obj.city,[obj.stations(iplus1).n,obj.stations(k).n]);
            new_links(4)=Path(obj.city,[obj.stations(jplus1).n,obj.stations(kplus1).n]);
            obj=new_concatenate([paths,new_links]);
            obj=obj.from_path_to_routes();

        end



        function obj=type2_insertion(obj,station,i,j,k,l)
            % type 2 insertion is described in the article by M.Gendreau,
            % A.Hertz,G.Laporte 'New insertion and post-optimization procedures
            % for the traveling salesman problem.'

            cond= k~=j && k~=j+1 && i~=j && l~=i && l~=i+1 && all(station.n~=[obj.stations.n]) ;
            if ~cond
                error("Type II insertion conditions not verified");
            end
            len=obj.n_stations-1;
            if obj.is_in_the_middle(i,j,k) 
                if obj.is_in_the_middle(j,i,l)
                obj=obj.invert();
                inds=mod(len-[i-1,j-1,k-1,l-1],len)+1;
                i=inds(1);
                j=inds(2);
                k=inds(3);
                l=inds(4);
                obj=type2_insertion(obj,station,i,j,k,l);
                return
                else
                    error("Type II insertion conditions not verified");
                end
            end

            kminus1=obj.plus(k,-1);
            lminus1=obj.plus(l,-1);
            iplus1=obj.plus(i,+1);
            jplus1=obj.plus(j,+1);

            
            paths=obj.break_path([i,j,kminus1,lminus1]);
            new_links=Path(obj.city,[obj.stations(i).n,station.n]);
            new_links(2)=Path(obj.city,[station.n ,obj.stations(j).n]);
            new_links(3)=Path(obj.city,[obj.stations(l).n,obj.stations(jplus1).n]);
            new_links(4)=Path(obj.city,[obj.stations(kminus1).n,obj.stations(lminus1).n]);
            new_links(5)=Path(obj.city,[obj.stations(iplus1).n,obj.stations(k).n]);            
            obj=new_concatenate([paths,new_links]);
            obj=obj.from_path_to_routes();

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %   unstringing methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=type1_unstringing(obj,i,j,k)
            % type 1 unstinging is described in the article by M.Gendreau,
            % A.Hertz,G.Laporte 'New insertion and post-optimization procedures
            % for the traveling salesman problem.'
            cond= i~=j  &&  obj.is_in_the_middle(i,j,k) && obj.plus(j,+1)~=i;%i~=k && j~=k;
            if ~cond
                error('type1_unstringing conditions not verfied')
            end
            len=obj.n_stations-1;

            iminus1=obj.plus(i,-1);
            breaking_points=[iminus1,i,k,j];
            paths=obj.break_path(breaking_points);

            iplus1=obj.plus(i,+1);
            kplus1=obj.plus(k,+1);
            jplus1=obj.plus(j,+1);

            new_links(3)=Path(); %preallocate
            new_links=Path(obj.city,[obj.stations(iminus1).n,obj.stations(k).n]);
            new_links(2)=Path(obj.city,[obj.stations(iplus1).n ,obj.stations(j).n]);
            new_links(3)=Path(obj.city,[obj.stations(kplus1).n,obj.stations(jplus1).n]);

            obj=new_concatenate([paths,new_links]);
            obj=obj.from_path_to_routes();

        end

        function obj=type2_unstringing(obj,i,j,k,l)
            % type 2 unstringing is described in the article by M.Gendreau,
            % A.Hertz,G.Laporte 'New insertion and post-optimization procedures
            % for the traveling salesman problem.'
            iminus1=obj.plus(i,-1);
            jminus1=obj.plus(j,-1);
            iplus1=obj.plus(i,+1);
            lplus1=obj.plus(l,+1);
            kplus1=obj.plus(k,+1);
            cond= i~=j && obj.is_in_the_middle(j,iminus1,k) && obj.is_in_the_middle(jminus1,k,l)...
                && lplus1~=i && kplus1~=i && l~=i && k~=i && jminus1~=i;
            if ~cond
                error('type2_unstringing conditions not verfied')
            end
            
            len=obj.n_stations-1;
            iminus1=obj.plus(i,-1);
            jminus1=obj.plus(j,-1);

            breaking_points=[iminus1,i,jminus1,l,k];
            paths=obj.break_path(breaking_points);

            iplus1=obj.plus(i,+1);
            lplus1=obj.plus(l,+1);
            kplus1=obj.plus(k,+1);

            new_links(4)=Path(); %preallocate
            new_links=Path(obj.city,[obj.stations(iminus1).n,obj.stations(k).n]);
            new_links(2)=Path(obj.city,[obj.stations(lplus1).n ,obj.stations(jminus1).n]);
            new_links(3)=Path(obj.city,[obj.stations(iplus1).n,obj.stations(j).n]);
            new_links(4)=Path(obj.city,[obj.stations(l).n,obj.stations(kplus1).n]);

            obj=new_concatenate([paths,new_links]);
            obj=obj.from_path_to_routes();


        end

        %
        % route optimization
        %

        function min_route=opt3(obj,max_it)
            % opt 3
            
            %
            % funzione d'appoggio   
            %
            function min_route=opt3_appo(obj)
                min_route=obj;
                %estraggo l valori da 2 a n-1 e divido la route in path
                n=1;
                breaks=zeros(nchoosek(obj.n_stations-1,3),3);
                for i=1:obj.n_stations-1
                    for j=i+1:obj.n_stations-1
                        for k=j+1:obj.n_stations-1
                            breaks(n,1)=i;
                            breaks(n,2)=j;
                            breaks(n,3)=k;
                            n=n+1;
                        end
                    end
                end

                %estraggo l valori da 2 a n-1 e divido la route in path

                for n=1:length(breaks)
                    paths= obj.break_path(breaks(n,:));

                    candidate_path_ordering= [1,2,3,4;1,3,2,4];

                    % this matrix tells if the conatenation is at the
                    % beginning or at the end of the path
                    directions=[zeros(2^2,1),...
                        [ones(2,1),[1;0];zeros(2,1),[1;0]],...
                        zeros(2^2,1)];

                    candidates= paths.concatenate(candidate_path_ordering,directions);
                    candidates=candidates.from_path_to_routes();

                    %min([candidates.cost])
                    if min([candidates.cost])< min_route.cost
                        costs=[candidates.cost];
                        k=find(costs==min(costs),1);
                        min_route=candidates(k);
                    end
                end
            end


            if nargin==1
                max_it=inf;
            end
            old_route=obj;
            k=1
            while k<max_it
                min_route=opt3_appo(old_route);
                if min_route.cost==old_route.cost
                    return
                end
                k=k+1
                old_route=min_route;
            end
        end



        %
        % support methods
        %

        function stations=nbhd(obj,station,p)
            %evaluate the neighborood of a station
            route_stations=[obj.stations.n];
            route_stations=route_stations(1:end-1);
            dist=obj.city.distances(station.n,route_stations);
            [~,inds]=sort(dist);
            len=length(inds);
            if dist(inds(1))==0 %significa che la stazione di cui vedere il vicinato è nella route
                stations=[obj.stations(inds(2:min([len,p+1])))];
                return
            end
            stations=[obj.stations(inds(1:min([len,p])))];
        end


        function station=closest_station(obj)
            %CLOSEST_STATION finds the closest station to a route in the sense of the
            %minimum of the minimum distance to the stations of the route
            
            dist=obj.city.distances;
            %boolean array of the stations in the city that are in the
            %route
            points_in_route=zeros(obj.city.n_stations,1);
            points_in_route([obj.stations.n])=1;
            minimum_distance=min(min(dist(find(points_in_route),find(~points_in_route))));
            [r_candidates,c_candidates]=ind2sub(size(dist)  ,find(dist==minimum_distance));

            ind=find(ismember(r_candidates,[obj.stations.n]) & ~ismember(c_candidates,[obj.stations.n]),1);
            ind=c_candidates(ind);
            station=obj.city.stations(ind);

        end


        function bool=is_in_the_middle(obj,i,j,k)
            %this function evaluate if k is in the middle of i,j
            bool= (i<j && k>i && k<j) || (j<i && (k>i || k<j));
        end

        function out=remove(obj,station)
            ind = find([obj.stations.n]==station.n,1);    %questo è l'indice sulla strada
            out=route();
            out.city=obj.city;
            ind_plus_1=min(obj.n_stations,ind+1);
            out.stations=obj.stations([1:ind-1,ind_plus_1:end]);
            out.demand=obj.demand-station.demand;
            out.cost= obj.cost +...
                obj.city.distances(obj.stations(obj.plus(ind,-1)).n,obj.stations(obj.plus(ind,1)).n)-...
                obj.city.distances(obj.stations(obj.plus(ind,-1)).n,station.n)-...
                obj.city.distances(station.n,obj.stations(obj.plus(ind,1)).n);

                
            if obj.cost<0
                error('errore')
            end

            out.n_stations = obj.n_stations-1;

        end

        function out=copy(obj,out)
            % copy is needed since route is an handle class
            if nargin<2
                out=route();
            end
            out.city=obj.city;
            out.stations=obj.stations;
            out.demand=obj.demand;
            out.cost= obj.cost;
            out.n_stations = obj.n_stations;

        end

        function out=embedd(obj,city)
            % if possible it embedds a route in a city into another city 
            ind=zeros(1,obj.n_stations);
            for i=1:obj.n_stations
                for v=city.stations
                    if v.x==obj.stations(i).x &&...
                           v.y==obj.stations(i).y &&...
                           v.demand==obj.stations(i).demand &&...
                           all(v.category==obj.stations(i).category)
                        
                        ind(i)=v.n;
                    end
                end
            end
            if any(ind==0)
                error('non si può fare embedding di città')
            end
            ind;
            out=route(city,ind);
        end

        function fig=plot(obj)
            % Plot: Visualization method
            fig=plot(obj(1).city);
            fig;
            hold on;
            for i=1:length(obj)
                x=[obj(i).stations.x]; 
                y=[obj(i).stations.y];
                a=plot(x(2:end-1),y(2:end-1));
                plot(x([1,2]),y([1,2]),'--',"Color",a.Color)
                plot(x([end-1,end]),y([end-1,end]),'--',"Color",a.Color)
            end
            
            hold off;
        end
    end
end