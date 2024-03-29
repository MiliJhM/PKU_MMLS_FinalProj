% filename main.m 
function T = main(pliststep, plistrang, plistifergo)
% Main function
% Grid traversal of Gene drive related parameters

%% Section 1:: Finished
% Initialize Paramlist for Cycle

% 本节所有结构体中参数仅与基因驱动过程有关，其余均采用环境拟合参数

plistDef = struct(...
    'drive_conversion', 0.8, ...
    'drive_fitness', 0.9, ...
    'release_rate', 0.001, ...
    'germline_resistance_forming', 0.1, ...
    'dd_mothertal_inheritance', 0.05,...
    'log_immunity_speed', 0,...
    'reducehtm', 0.5 ...
    );

% Cycle:
% 对前面的所有参数相关设置，对每个参数转化为参数的可行取值列表
for field = fieldnames(plistrang)'
    field = field{1};
    if plistifergo.(field)
        min_val = plistrang.(field)(1);
        max_val = plistrang.(field)(2);
        step = pliststep.([field '_step']);
        plist.(field) = min_val:step:max_val;
    else
        plist.(field) = plistDef.(field);
    end
    
end

%% Section 2 :: Finished
% Const defination
% Struct phumanlist/pmosquitolist: 
% List of parameters obtained by fitting real data
phumanlist = struct( ...
    'human_recovery', {0.030673973459249}, ...
    'mtoh', {0.110070995896335}, ...
    'immunity_gain_rate', {0.086442942066365}, ...
    'immunity_losing_rate', {0.001545980207871}, ...
    'b1', {0.801405190261370},...
    'shape',{2.013061848936364}...
    );

pmosquitolist = struct( ...
    'capacity', {29}, ...
    'reproduction_rate', {0.5}, ...
    'htom', {0.2} ...
    );

% Real data used for fitting
EIRS=[3.748022899,4.209665958,20.302014,19.46226874,37.85310084,68.37694864,41.62689645,67.65891383,92.87530606,128.8428103,126.1500306,157.4671712,198.6464579,209.4149771,240.2292953,169.5482093,159.1395626,162.5352447];
EIRS=reshape(EIRS,length(EIRS),1);
prevalences=[0.02769617,0.167517702,0.098683032,0.283678023,0.25356346,0.268620741,0.612798931,0.511697506,0.470825206,0.685937987,0.617101703,0.59559107,0.580532176,0.552568192,0.602044422,0.806399466,0.892445223,0.944072032];
prevalences=reshape(prevalences,length(prevalences),1);

% Const defination end


%% Section 3
% Main Cycle
% Perform parallel operation on the previously obtained Gene drive parameter list to obtain the average Prev minimum value for each parameter
paramMatrix = {plist.drive_conversion, plist.drive_fitness, plist.release_rate, plist.germline_resistance_forming, plist.dd_mothertal_inheritance, plist.log_immunity_speed, plist.reducehtm};
param_combinations = combvec(paramMatrix{:});
results = zeros(length(param_combinations), 1);

parfor i = 1:length(param_combinations)
    results(i) = spacial('Results/findings/param',make_paramlist(param_combinations(:,i)),phumanlist,pmosquitolist,reshape(1:250,250,1),reshape(1:250,250,1),true,false,false,false);
end

%% Section 4
% Data collection
% Collect the previously obtained data for plotting and analysis
driveConversion = param_combinations(1,:)' ;
driveFitness = param_combinations(2,:)';
releaseRate = param_combinations(3,:)';
germlineResistanceForming = param_combinations(4,:)';
mothertalInheritance = param_combinations(5,:)';
logImmunitySpeed = param_combinations(6,:)';
reduceHtm = param_combinations(7,:)';

T = table(driveConversion,driveFitness,releaseRate,germlineResistanceForming,mothertalInheritance,logImmunitySpeed, reduceHtm, results);
dir = sprintf('Results/result%d%d%d%d%d%d%d.csv', plistifergo.drive_conversion, plistifergo.drive_fitness, plistifergo.release_rate, plistifergo.germline_resistance_forming, plistifergo.dd_mothertal_inheritance, plistifergo.log_immunity_speed, plistifergo.reducehtm);
writetable(T, dir)

end