


function [S_min,S_tilde_min,F1_min,alpha,fv] = SEARCH(W,q,p1,p2,theta_min,theta_max,g,h,n_max...
    ,vrp,alpha,fv)
%SEARCH is a tabu search algorithm for vrp
%   input:
%       W = non empty susbset of stations that are allowed to be moved
%       q = number of vertices of W that are candidate for reinsertion into
%       another route
%       p1 = number of neighbors to be allowed to be inserted
%       p2 = nbhds size used in GENI
%       theta_min,theta_max = number iteration of tabu memory
%       g = scaling factor
%       h = frequency of alpha and beta updates
%       n_max = max number of non imprevement steps 




%
% step 0 (initialization)
%

S=vrp.solution;
m_tilde=S.n_vehicles; 
m=vrp.n_vehicles;
city = vrp.city;
dist=city.distances;
Q=vrp.capacity;
S_min=VrpSolution.empty;
if S.is_feasible(Q)
    S_min=S;
end
S_tilde_min=VrpSolution.empty;

Q=vrp.capacity;
t=1;
tabuMoves=Tabu.empty;

F1= @(sol) sol.cost;
F2= @(sol,alpha) F1(sol) +...
    alpha * sum ( max([sol.routes.demand]-Q,0) ); 

F1_min=inf; %F1(S);
F2_min=inf; %F2(S,alpha);

non_decreased=0;
usedUS=0;
lambda_max=0;
F_2=F2(S,alpha);


R_minus_v(q)=route();
last_h_solutions_feasibility=zeros(1,h);



while(true)

    init=tic;

    %
    % step 1 (vertex selection)
    %

    selected_stations=randsample(W,q,false);

    parfor v=1:length(selected_stations)

        %consider all potential moves of v from its current route into another
        %route containing no city if possible or at least one of the p1 nearest
        %neighbors


        Rr=find_route(S,selected_stations(v)); % trova la route in cui sta una soluzion
        p2=max(p1,floor(Rr.n_stations/2));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Genero le routes candidate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        candidate_routes=route.empty;
        if m_tilde<m
            candidate_routes = route(vrp.city,[1,1]);
        end

        %calcolo p1 neigbors
        [~,v_nbhd]=sort(dist(selected_stations(v).n,:));
        index_depot=find(v_nbhd==1);
        v_nbhd=v_nbhd([1:index_depot-1,index_depot+1:end]);        
        v_nbhd=v_nbhd(2:p2+1); % nel nbhd non ci deve essere il depot
        
        candidate_routes_ind=zeros(1,length(S.routes));
        for r=1:length(S.routes) %routes è l'attributo con l'array di routes
            for i=v_nbhd
                if ~isequal(S.routes(r),Rr) && any(i==[S.routes(r).stations.n])
                    candidate_routes_ind(r)=1;
                end
            end
        end

        candidate_routes=[candidate_routes,S.routes(logical(candidate_routes_ind))];


        %
        % (a)   remove v from Rr and compute its insertion cost into Rs using
        %       GENI with parameter p2 end determines S'
        %

        R_minus_v(v)=Rr.remove(selected_stations(v));
        
        if isempty(candidate_routes)
            %if there are no candidate routes extract new stations and
            % repeat
            Frmin(v)=inf;
            continue
        end


        S_post=VrpSolution.empty; %elimino interazione precedente
        S_post(length(candidate_routes))=VrpSolution(); %inizializzo
        F=zeros(length(candidate_routes),1);

        for r=1:length(candidate_routes)
            old_routes=S.routes(S.routes~=Rr & S.routes~=candidate_routes(r));
            Rs=copy(candidate_routes(r));
            Rs=GENI_INSERTION(selected_stations(v),Rs,p2);
            S_post(r)=VrpSolution([old_routes,Rs,R_minus_v(v)]);


            %
            % (b)   if the move is tabu, it is disregarded unless:
            %           - S' is feasible and F1(S')<F1_min
            %           - S' infeeasible and F2(s')<F2_min


            if tabuMoves.is_tabu(selected_stations(v),candidate_routes(r)) && ...
                    ~((S_post(r).is_feasible(vrp.capacity) && F1(S_post(r))<F1_min) ||...
                    (~S_post(r).is_feasible(vrp.capacity) && F2(S_post(r),alpha)<F2_min))
                
                %discard solution if the move is tabu
                F(r)=inf;
           
            else
                %
                % (c)   otherwise S' is assigned a value F(S')=F2(S') se F2(S')<F2(S)
                %       o F(S')=F2(S')+lambda_max \sqrt(m) g f_v
                %           -lambda_max =  is the largest observed absolute difference
                %           between the values of F2(S)obtained a two successive
                %           interation
                %           - f_v number of times vertex v has been moved divided by t
                %

                if F2(S_post(r),alpha)<F2(S,alpha)
                    F(r)=F2(S_post(r),alpha);
                else
                    F(r)=F2(S_post(r),alpha) + lambda_max *m^(1/2) *g * fv(selected_stations(v).n);
                end
            end

        end

        [minimum,ind]=min(F);
        Frminroutes(v)=candidate_routes(ind); %ho il puntatore di dove sono arrivato
        old_routes_min{v}=S.routes(S.routes~=Rr & S.routes~=candidate_routes(ind)); %lista delle routes non cambiate nel minimo
        Frmin(v)=minimum;
        S_post_v(v)=S_post(ind);


    end

    %
    % step 3 (Identification of best move)
    %
    % The candidate move yelding the least value of F and solution S_bar is
    % identified
    %
    [~,ind]=min(Frmin);

    if Frmin(ind)==inf %se non trovo niente di possibile ri estraggo i v
        warning('non ho trovato candidate routes')
        non_decreased=non_decreased+1;
        if non_decreased>n_max
            return
        end

        t=t+1;
        continue
    end

    S_bar=S_post_v(ind);
    changed_station=selected_stations(ind);
    Rr=find_route(S,changed_station); %puntatore della route da cui me ne sono andato (prima)
    Rr_new=R_minus_v(ind); %puntato della route da cui me ne sono andatao (dopo)
    Rs=Frminroutes(ind); %puntatore della route in cui sono arrivato (prima)
    Rs_new=find_route(S_bar,changed_station); %puntatore della route in cui sono arrivato (dopo)
    R_old=old_routes_min{ind}; %puntatori delle routes non cambiate


    %
    % Step 4 (next solution)
    %
    % The move identified in Step 3 is not necessarily implemented. It may
    % indeed be advantageus to attempt to improve S by applying to each
    % individual route of D the US post-optimizaion procedure. Solution S=S_bar
    % unless the three condition are satisfied:
    % (a) F2(S_bar)>F_2(S)
    % (b) S is feasible
    % (c) US has not been used at iteration t-1
    % In such case S is obtained by applying US
    %

    %improve S
    if usedUS==0
        route_appo=route.empty;
        route_appo(length(S.routes))=route();
        ind=mod(t-1,h)+1;
        if ind==h
            parfor r=1:length(S.routes)

                %S.routes(r)=
                route_appo(r)=US(S.routes(r));
                %route_appo(r)=S.routes(r);
            end
            S.routes=route_appo;
        end
    end


    %%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %S_OLD=S;
    %for k=1:length(S.routes)
    %    rou(k)=S.routes(k).copy()
    %end
    %S_OLD.routes=rou;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    if F2(S_bar,alpha)>F2(S,alpha) &&...
            S.is_feasible(Q) &&...
            usedUS==0
        usedUS=1;
    else
        usedUS=0;
        %voglio però anche che non siano cambiati i puntatori
        Rr=copy(Rr_new,Rr);
        Rs=copy(Rs_new,Rs);
        S=VrpSolution([R_old,Rr,Rs]);



    end


    %
    % Step 5 (Update)
    %
    % If the US procdure has not been used in Step 4 and vertex v has been
    % moved from Rr to Rs, reinsertinf v in Rr is considered tabu until t+theta
    % where theta is integer Unif~[Theta_min, Theta_max]
    % t=t+1
    % Update F1_min, F2_min, S_min, S_tilde_min, lambda_max, m , fv

    if ~usedUS
        tabuMoves(end+1)=Tabu(changed_station,Rr,theta_min,theta_max);
    end

    t=t+1;
    tabuMoves=tabuMoves.next();

    S.routes=S.routes([S.routes.n_stations]~=2); %elimino routes inutili
    for k=1:length(S.routes)
        if S.routes(k).stations(1).n~=1
            ind=find([S.routes(k).stations.n]==1,1);
            S.routes(k).stations=S.routes(k).stations([ind:S.routes(k).n_stations,2:ind]);
        end
    end

    non_decreased=non_decreased+1;
    if F1(S)<F1_min && S.is_feasible(Q)
        F1_min=F1(S);
        S_min=copy(S);
        non_decreased=0;
    end

    if F2(S,alpha)<F2_min
        F2_min=F2(S,alpha);
        S_tilde_min=copy(S);
        non_decreased=0;
    end

    lambda_max=max(lambda_max,abs(F_2-F2(S,alpha)));
    F_2=F2(S,alpha);

    m_tilde=S.n_vehicles;
    fv=(t-2)*fv; %numero di volte che ho spostato le stazioni
    fv(changed_station.n)=fv(changed_station.n)+1; %aggiungo uno alla stazione che ho spostato
    fv=fv/(t-1);


    %
    % Step 6 (Penalty adjustment)
    %
    % If t is a multiple of h adjust alpha and beta
    % Devo salvare le ultime h soluzioni

    ind=mod(t-1,h)+1;
    last_h_solutions_feasibility(ind)=S.is_feasible(Q);
    if ind==h
        if all(last_h_solutions_feasibility)
            alpha=alpha/2;
        elseif all(~last_h_solutions_feasibility)
            alpha=alpha*2;
        end
    end


    time=toc(init);
    fprintf('\n%d° iteration  \n%d° non decreasing iteration \n%f seconds for the iteration \n',t-1,non_decreased,time);

    station_in_sol=zeros(1,city.n_stations);
    for r=S.routes
        for station=r.stations
            station_in_sol(station.n)=station_in_sol(station.n)+1;
        end
    end


    %%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %if any(station_in_sol==0)
    %    error("c'è qualcosa che non va")
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %
    % Step 7 (termination check)
    % If F1_min e F2_min have not decreased for the last n_max iteration stop.
    % Otherwise go to step 1

    if non_decreased>n_max
        return
    end

end
end

