function sol = GENI_INSERTION(station,tour,p)
%GENI_INSERTION is a simplified version of the geni algorithm

city=station.city;


if  length(tour.stations)<=2
    sol=route(city,[1,station.n,1]);
    return
end

if  length(tour.stations)<=3
    tour_stations_inds=[tour.stations.n]
    sol=route(city,[tour_stations_inds(1:end-1),station.n,1]);
    return
end

%
%   step 1
%
%initialize tour

%new_stations=randsample(2:city.n_stations,2,'false');
%sol=route(city,[1,new_stations,1]);

points_in_tour=zeros(1,city.n_stations);
points_in_tour([tour.stations.n])=1;

if nargin==2
    p=3;
end
p_nb=zeros(city.n_stations);
for i=1:city.n_stations
    p_nb(i,[tour.nbhd(city.stations(i),p).n])=1;
end


%
%   step2
%


new_station=station;
sol=tour;


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

%   type 2 insert

for i=find(p_nb(new_station.n,:))   %i è un indice nella città
    i=find([sol.stations.n]==i,1);
    for j=find(p_nb(new_station.n,:))   %j è un indice nella città
        j=find([sol.stations.n]==j,1);

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


ns=[sol.stations.n];
ind=find([sol.stations.n]==1,1);
sol=route(sol.city,ns([ind:sol.n_stations,1:ind]));
end

