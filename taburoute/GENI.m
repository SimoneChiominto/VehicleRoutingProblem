function sol = GENI(stations,p)
% GENI is a tsp heuristic
%   input:
%       stations: stations we want to include in the cycle
%       p: numeber of neighbor we want to include
%   output: 
%       sol: object of class route that is an approximated solution of a
%            tsp 

original_city=stations(1).city;
city=City(stations);

if  length(stations)==2
    sol=route(city,[1,2,1]);
    return
end


%
%   step 1
%

%initialize tour
new_stations=randsample(2:city.n_stations,2,'false');
sol=route(city,[1,new_stations,1]);

points_in_tour=zeros(1,city.n_stations);
points_in_tour([1,new_stations])=1;

if nargin==1
    p=3;
end
p_nb=zeros(city.n_stations);
for i=1:city.n_stations
    p_nb(i,[sol.nbhd(city.stations(i),p).n])=1;
end



%
%   step2
%

for n=3:city.n_stations-1

    new_station=randsample(find(~points_in_tour),1);
    if length(find(~points_in_tour))==1
        new_station=find(~points_in_tour);
    end
    new_station=city.stations(new_station);


    candidates=route.empty;
    number_of_candidates=1;
    
    % type1 insert
    
    for i=find(p_nb(new_station.n,:))   %i è un indice nella città
        i=find([sol.stations.n]==i,1);
        for j=find(p_nb(new_station.n,:))   %j è un indice nella città
            j=find([sol.stations.n]==j,1);
            for next=[-1,1]   %due possibili ordinamenti della route
                ind=sol.plus(i,next);             %passo al successivo
                ind=sol.stations(ind).n;            %ritorno a indice su città
                for k=find(p_nb(ind,:))
                    k=find([sol.stations.n]==k,1);
                %for k=1:n
                    try
                        candidates(number_of_candidates)=sol.type1_insertion(new_station,i,j,k);
                        number_of_candidates=number_of_candidates+1;
                    catch
                        continue;
                        
                    end
                end
            end
        end
    end

    % type2 insert

    for i=find(p_nb(new_station.n,:))   %i è un indice nella città
        i=find([sol.stations.n]==i,1);
        for j=find(p_nb(new_station.n,:))   %j è un indice nella città
            j=find([sol.stations.n]==j,1);
            % v_ki N v_i+1              
            for next=[-1,1]   %due possibili ordinamenti della route                
                ind=sol.plus(i,next);             %passo al successivo
                ind=sol.stations(ind).n;            %ritorno a indice su città
                for k=find(p_nb(ind,:))
                    k=find([sol.stations.n]==k,1);                    
                    ind=sol.plus(j,next);             %passo al successivo
                    ind=sol.stations(ind).n;            %ritorno a indice su città
                    for l=find(p_nb(ind,:))
                        l=find([sol.stations.n]==l,1);
                        try
                            candidates(number_of_candidates)=sol.type2_insertion(new_station,i,j,k,l);
                            number_of_candidates=number_of_candidates+1;
                        catch
                            continue;
                        end
                    end
                end
            end
        end
    end

    [~,I]=min([candidates.cost]);
    sol=candidates(I);

    points_in_tour([sol.stations.n])=1;
    %update p_nbhds
    p_nb=zeros(city.n_stations);
    for i=1:city.n_stations
        p_nb(i,[sol.nbhd(city.stations(i),p).n])=1;
    end

    
end

ns=[sol.stations.n];
ind=find([sol.stations.n]==1,1);
sol=route(sol.city,ns([ind:sol.n_stations,1:ind]));
sol=embedd(sol,original_city);

end
