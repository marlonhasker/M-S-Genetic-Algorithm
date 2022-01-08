clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in Data:                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix containing detection chances:
S = readmatrix('S.txt');
% Matrix containing travel times (in minutes):
T = readmatrix('T.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These parameters are fixed:      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the survival chance is less than this, the uav will be considered lost
cuttoff_Survival = 0.3;
% If travel time is larger than this, the target will have moved:
cuttoff_Time = 24;

n = 20; % number of rows
m = 10; % number of columns
t = m * n + 2; % node number of destination

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Genetics Algorithm:              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
population_size = 1000;
population = cell(population_size, 3);
for i = 1:population_size
    population(i,:) = generate_random_route(n, m, t, T, S, cuttoff_Survival, cuttoff_Time);
end

% Run generations
num_generations = 200;
avg_fitness = zeros(num_generations,3);
death_penalty = 0.15; % Percentage of least fit population to be removed.
copy_best = 0.10; % Percentage of fittest population to be directly copied.
mutation_chance = 0.12; % Chance of mutation per column
for i=1:num_generations
    population = sortrows(population,[2 3],{'descend', 'ascend'}); % Prioritize survival, then time
    avg_fitness(i,1) = population{1,2}; % Store best fitness
    avg_fitness(i,2) = mean([population{1:floor(population_size/10),2}]); % Top 10% fitness
    avg_fitness(i,3) = mean([population{:,2}]); % Average fitness

    population = population(1:floor(population_size*(1-death_penalty)),:); % Apply Death Penalty.
    if mod(length(population),2) % Make sure population size is even for crossover.
        population(end,:) = [];
    end
    upper_class = population(1:ceil(population_size*copy_best),:);

    parent_pairs = reshape(randperm(length(population)),[],2); % Make pairs of parents
    child_population = cell(length(population),3);
    for k=1:length(parent_pairs) % Make children
        pair = [population(parent_pairs(k,1),1) population(parent_pairs(k,2),1)];
        child_population(k*2-1:k*2,:) = generate_children(pair, mutation_chance, n, m, T, S, cuttoff_Survival, cuttoff_Time);
    end
    child_population = sortrows(child_population,[2 3],{'descend', 'ascend'});
    child_population(end-length(upper_class)+1:end,:) = []; % Remove least fit children

    num_immigrants = population_size-length(upper_class)-length(child_population);
    immigrant_population = cell(num_immigrants,3);
    for k=1:num_immigrants % Generate immigrants to restore population size
        immigrant_population(k,:) = generate_random_route(n, m, t, T, S, cuttoff_Survival, cuttoff_Time);
    end
    population = [upper_class; child_population; immigrant_population];
end

population = sortrows(population,[2 3],{'descend', 'ascend'}); %
best_route = population(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots:                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,2,1)
hold all
plot(avg_fitness(:,1));
plot(avg_fitness(:,2));
plot(avg_fitness(:,3));
title('Fitness per generation')
xlabel('Generation')
ylabel('Chance of survival')
ylim([0.3 1.0])
legend({'Best route','Top 10%','Entire population'},'Location','southeast')
grid on
hold off 


% The following 8 lines ensure a nice structured graph of T & S:
horizontalspacing = 2;
x = zeros(1,t);
for i=1:m
    x((i-1)*n+2:i*n+1) = i*horizontalspacing;
end
x(end) = (m+1)*horizontalspacing;
y = [0 repmat(linspace(-2.5, 2.5, n), 1, m) 0];
z = zeros(1,length(x));

subplot(1,2,2)
g1 = graph(T, 'upper','omitselfloops');
hold all
p = plot(g1,'XData',x,'YData',y,'ZData',z);
view(2)
highlight(p,best_route{1}, 'LineWidth', 5, 'EdgeColor', 'g')
title('Optimal Route');
set(gca,'XColor', 'none','YColor','none')
set(gcf,'units','inch','position',[4,4,20,8])
route_description = "Node" + sprintf(' %d →',best_route{1}) + ...
    sprintf('\nChance of survival: %.1f%%, Time: %d', population{1,2}*100, population{1,3});
text(0,-3,route_description)
hold off  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions:                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function route = generate_random_route(n, m, t, T, S, cuttoff_Survival, cuttoff_Time)
    path = [1 zeros(1,m) t];
    viable = false;
    while ~viable
        for i=2:m+1
            path(i) = 1 + (i-2)*n + randi(n);
        end
        f = fitness(path, T, S);
        viable = will_survive(f, cuttoff_Survival, cuttoff_Time);
    end
    route = {path, f(1), f(2)};
end

function candidate = mutate_route(route, mutation_chance, n, m, T, S, cuttoff_Survival, cuttoff_Time)
    viable = false;
    while ~viable
        candidate = route;
        for i=2:m+1
            if rand <= mutation_chance
                candidate{1}(i) = 1 + (i-2)*n + randi(n);
            end
        end
        f = fitness(candidate{1}, T, S);
        viable = will_survive(f, cuttoff_Survival, cuttoff_Time);
    end
    candidate = {candidate{1}, f(1), f(2)};
end

function children = generate_children(parents, mutation_chance, n, m, T, S, cuttoff_Survival, cuttoff_Time)
    num_children = 2;
    children = cell(num_children, 3);
    childpath = cell(3);
    for i=1:num_children
        viable = false;
        while ~viable
            a = randi(2)-1;
            p1 = a+1;
            p2 = ~a+1;
            splicing_point = 1 + randi(m);
            childpath = [parents{p1}(1:splicing_point) parents{p2}(splicing_point+1:end)];
            f = fitness(childpath, T, S);
            viable = will_survive(f, cuttoff_Survival, cuttoff_Time);
        end
        child = {childpath, f(1), f(2)};
        children(i,:) = mutate_route(child, mutation_chance, n, m, T, S, cuttoff_Survival, cuttoff_Time);
    end
end

function fitness = fitness(path, T, S)
    time = 0;
    chance = 1.0;
    for i=1:length(path)-1
        time = time + T(path(i),path(i+1));
        chance = chance * (1 - S(path(i),path(i+1)));
    end
    fitness = [chance, time];
end

function viable = will_survive(fitness, cuttoff_Survival, cuttoff_Time)
    viable = true;
    if fitness(1) <= cuttoff_Survival || fitness(2) > cuttoff_Time
        viable = false;
    end
end




