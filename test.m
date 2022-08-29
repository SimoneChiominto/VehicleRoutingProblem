clear 
addpath('classes')
addpath('taburoute')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% insert input
file='VRP_instances/E030-03g.dat';
n_vehicles=3;
capacity=4500;
test_file='TEST_results/TEST_E030-03.vrp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

city=City();
city=city.initialize_from_file(file);

vrp=VRP(city,n_vehicles,capacity);
[vrp_CFRS,time_CFRS]=vrp.solve('CFRS');
[vrp_CW,time_CW]=vrp.solve('Clarke-Wright');
[vrp_tabu,time_tabu]=vrp.solve('tabu-search');


%save(test_file)


