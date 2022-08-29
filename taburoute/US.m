function sol = US(sol,p)
%US is an optimization algorithm for tsp solutions
%   input:
%       sol: tsp solution to be optimized
%       p: number of neighbors to be considered
%   output:
%       sol: tsp solution after US optimization

if sol.n_stations<=4
    return
end


%
%   Step 1
%
tau_min=sol;
z_min=sol.cost;
t=1;

%
%   Step 2
%

if nargin==1
    p=3;
end
p_nb=zeros(sol.city.n_stations);
for i=1:sol.city.n_stations
    p_nb(i,[sol.nbhd(sol.city.stations(i),p).n])=1;
end

while t<sol.n_stations

    t_station=sol.stations(t);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UNSTRINGING mi dà una route levato un nodo
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %array di soluzioni senza il nodo t
    sol_sans_t=unstringing(tau_min,t,p_nb); %array di soluzioni senza il nodo t

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STRINGING mi connette la route
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sol_candidates=stringing(sol_sans_t,t_station,p); %array di candidate soluzioni


    tau_appo=sol_candidates(find(min([sol_candidates.cost])==[sol_candidates.cost],1));
    z_appo=min([sol_candidates.cost]);

    if z_appo<z_min
        tau_min=tau_appo;
        z_min=z_appo;
        t=1;
    else
        t=t+1;
    end
end

sol=tau_min;

end



function candidates=unstringing(sol,i,p_nb)

candidates=route.empty;
number_of_candidates=1;

for next=[-1,1]   %due possibili ordinamenti della route
    ind=sol.plus(i,next);             %passo al successivo
    ind=sol.stations(ind).n;            %ritorno a indice su città
    for j=find(p_nb(ind,:))   %j è un indice nella città
        j=find([sol.stations.n]==j,1);
        ind=sol.plus(i,-next);             %passo al successivo
        ind=sol.stations(ind).n;            %ritorno a indice su città
        for k=find(p_nb(ind,:))
            k=find([sol.stations.n]==k,1);
            try
                candidates(number_of_candidates)=sol.type1_unstringing(i,j,k);
                number_of_candidates=number_of_candidates+1;
            catch
                continue;
            end
        end
    end


    ind=sol.plus(i,next);             %passo al successivo
    ind=sol.stations(ind).n;            %ritorno a indice su città
    for j=find(p_nb(ind,:))   %j è un indice nella città
        j=find([sol.stations.n]==j,1);
        ind=sol.plus(i,-next);             %passo al successivo
        ind=sol.stations(ind).n;            %ritorno a indice su città
        for k=find(p_nb(ind,:))
            k=find([sol.stations.n]==k,1);
            ind=sol.plus(k,next);             %passo al successivo
            ind=sol.stations(ind).n;            %ritorno a indice su città
            for l=find(p_nb(ind,:))
                l=find([sol.stations.n]==l,1);
                try
                    candidates(number_of_candidates)=sol.type2_unstringing(i,j,k,l);
                    number_of_candidates=number_of_candidates+1;
                catch
                    continue;
                end
            end
        end
    end
end


end



function candidates=stringing(sols,new_station,p)

sol=sols(1);

% i nhds non vanno più bene e devono essere ricalcolati
p_nb=zeros(sol.city.n_stations);
for i=1:sol.city.n_stations
    p_nb(i,[sol.nbhd(sol.city.stations(i),p).n])=1;
end

%voglio che ci siano candidati per ogni soluzione proposta
candidates=route.empty;
number_of_candidates=1;

for sol=sols

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
end



end




