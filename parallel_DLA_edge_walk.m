function [radiuses, max_probs, edge_arr_sizes] = ...
    DLA_edge_walk_par(N, edge_walk_end_prob, min_num_of_simulated_particles,...
    max_iterations, draw, cut_zeros)

MIN_RADIUS = parallel.pool.Constant(10);
DELTA_R = parallel.pool.Constant(2);
PARTICLE_SYMB = parallel.pool.Constant(2);

if draw
cc = [1 1 1; 0 0 1]; % Blue and yellow (RGB values)
colormap(cc); % Put these two colors into the color map
set(0,'defaultaxesfontsize',20) % Make the default font size bigger
end

grid = zeros(N,N); % Domain where growth can occur.
% Empty sites are 0, cluster sites 2, placeholders are 3
grid(N/2,N/2) = 2; % Seed of the cluster
edge_arr = initial_edge(N);
rmax = 1; % Outer edge of the cluster
rad = 2; % particles are released from radius rmax + rad;

if draw
h = image(grid);
axis image
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
end

% Keep the cluster growing until it touches any edge of the
% simulation cell.
iteration = 0;
radiuses = zeros(1, max_iterations);
max_probs = zeros(1, max_iterations);
edge_arr_sizes = zeros(1, max_iterations);

% while the cluster isn't big enough or not enough particles simulated
while (all(grid(:,1)==0) && all(grid(:,N)==0) && ...
        all(grid(1,:)==0) && all(grid(N,:)==0)) || iteration > max_iterations
    iteration = iteration +1;
    radiuses(iteration) = rmax;
    edge_arr_sizes(iteration) = numel(edge_arr);
    r = max(rmax+rad, 10); % Pick a distance just beyond the edge of the cluster 
    if r>N/2; break; end % Stop when the cluster gets too big

    num_simulated_particles = max(min_num_of_simulated_particles, ...
       3 * edge_arr_sizes(iteration)^2); 
    chosen_edges = struct('x', [], 'y', [], 'from_x', [], 'from_y', [], 'valid', [], 'counter', []);
    
    collision_points = struct('x', [], 'y', [], 'from_x', [], 'from_y', []);
    
    %run all the particles in parallel until they hit the aggregate
    parfor(i=1:num_simulated_particles)
        [x, y, fr_x, fr_y] = simulate_single_particle(grid, r);
        collision_points(i).x = x;
        collision_points(i).y = y;
        collision_points(i).from_x = fr_x;
        collision_points(i).from_y = fr_y;
    end
    % after all the particles stuck, let the walk along the edge and
    % see where they end up
    for i=1:num_simulated_particles
        x = collision_points(i).x;
        y = collision_points(i).y;
        fr_x = collision_points(i).from_x;
        fr_y = collision_points(i).from_y;
        [new_edge, edge_arr] = ...
            choose_new_edge(edge_arr, x ,y ,fr_x, fr_y, edge_walk_end_prob);
        chosen_edges(i) = new_edge;
    end
    % update our output and choose an edge to add
    max_probs(iteration) = max(arrayfun(@(k) edge_arr(k).counter,...
                1:numel(edge_arr)))/num_simulated_particles;
    for j=1:numel(edge_arr)
        edge_arr(j).counter = 0;
    end
    index = randi(num_simulated_particles);
    edge_to_add = chosen_edges(index);
    
    % update the edge array after selecting a single edge to add
    [grid, edge_arr] = insert_edge(grid, edge_arr, edge_to_add, N);
    if ((edge_to_add.x-N/2)^2+(edge_to_add.y-N/2)^2)>rmax^2
        rmax = sqrt((edge_to_add.x-N/2)^2+(edge_to_add.y-N/2)^2);
    end

    if draw
    set(h,'cdata',grid)
    drawnow
    end
end
if cut_zeros
radiuses = radiuses(1,1:iteration);
max_probs = max_probs(1,1:iteration);
edge_arr_sizes = edge_arr_sizes(1:iteration);
end
end


function [grid, edge_arr] = insert_point(grid, edge_arr, x, y, fr_x, fr_y, prob)
[new_edge, ~] = choose_new_edge(edge_arr, x, y, fr_x, fr_y, prob);

new_x = new_edge.x; new_y = new_edge.y;
if grid(new_x, new_y) == 2
    error = 'dk'
end
grid(new_x, new_y) = 2;

edge_arr = update_edge_and_erase_holes(edge_arr, new_edge, grid);
end

function [grid, edge_arr] = insert_edge(grid, edge_arr, new_edge, N)
new_x = new_edge.x; new_y = new_edge.y;
if grid(new_x, new_y) == 2
    error = 'dk'
end
grid(new_x, new_y) = 2;

edge_arr = update_edge_and_erase_holes(edge_arr, new_edge, grid, N);
end

function [new_edge, edge_arr] = choose_new_edge(edge_arr, x, y, fr_x, fr_y, prob)
radius_sq = geornd(prob);
radius = floor(sqrt(radius_sq));
step_size = randi([-radius, radius]);

edge_point_eq = @(i) edge_arr(i).x == x && edge_arr(i).y == y ...
    && edge_arr(i).from_x == fr_x && edge_arr(i).from_y == fr_y ;
tf = arrayfun(edge_point_eq, 1:numel(edge_arr));
index = find(tf);
new_edge_index = mod(index - 1 + step_size, numel(edge_arr)) + 1;
new_edge = edge_arr(new_edge_index);
while ~new_edge.valid
    if step_size > 0 % when step_size=0 edge is always valid
        step_size = step_size -1;
    else
        step_size = step_size +1;
    end
    new_edge_index = mod(index - 1 + step_size, numel(edge_arr)) + 1;
    new_edge = edge_arr(new_edge_index);
end
edge_arr(new_edge_index).counter = edge_arr(new_edge_index).counter +1;
end

function [edge_arr] = update_edge(edge_arr, new_edge, grid)
x = new_edge.x; y=new_edge.y; fr_x = new_edge.from_x; fr_y = new_edge.from_y;
%first add new edges
edge_point_eq = @(i) edge_arr(i).x == x && edge_arr(i).y == y ...
    && edge_arr(i).from_x == fr_x && edge_arr(i).from_y == fr_y ;
tf = arrayfun(edge_point_eq, 1:numel(edge_arr));
index = find(tf);
edge_arr = add_edges(edge_arr, index, grid);

% erase/mark as not valid edges that point to the newly populated square
x = new_edge.x; y = new_edge.y;
to_erase_fn = @(i) edge_arr(i).x == x && edge_arr(i).y == y;
tf = arrayfun(to_erase_fn, 1:numel(edge_arr));
to_erase_indices = find(tf);

to_erase_tf = logical(zeros(1, numel(to_erase_indices)));
for i=1:numel(to_erase_indices)
    current = to_erase_indices(i);
    next = mod(current, numel(edge_arr)) +1;
    prev = mod(current -2, numel(edge_arr)) +1;
    px = edge_arr(prev).x; py = edge_arr(prev).y;
    nx = edge_arr(next).x; ny = edge_arr(next).y;
    if max(abs(px-nx), abs(py-ny)) >= 2  && false
        edge_arr(current).valid = false;
    else
        to_erase_tf(i) = true;
    end
end
to_erase_indices = to_erase_indices(to_erase_tf);
edge_arr(to_erase_indices) = [];
% edge_arr(min(to_erase_indices):max(to_erase_indices)) = [];
end

function [edge_arr] = update_edge_and_erase_holes(edge_arr, new_edge, grid, N)
x = new_edge.x; y=new_edge.y; fr_x = new_edge.from_x; fr_y = new_edge.from_y;
% first find the previous edge before modifying
edge_point_eq = @(i) edge_arr(i).x == x && edge_arr(i).y == y ...
    && edge_arr(i).from_x == fr_x && edge_arr(i).from_y == fr_y ;
tf = arrayfun(edge_point_eq, 1:numel(edge_arr));
index = find(tf);
prev_index = mod(index-2, numel(edge_arr))+1;
prev_edge = edge_arr(prev_index);
while prev_edge.x == x && prev_edge.y == y
    prev_index = mod(prev_index-2, numel(edge_arr))+1;
    prev_edge = edge_arr(prev_index);
end

% size_before = numel(edge_arr);

% erase
to_erase_fn = @(i) edge_arr(i).x == x && edge_arr(i).y == y;
tf = arrayfun(to_erase_fn, 1:numel(edge_arr));
to_erase_indices = find(tf);
if numel(to_erase_indices) == 1 %prev_index < min(to_erase_indices) || prev_index > max(to_erase_indices)
    edge_arr(min(to_erase_indices):max(to_erase_indices)) = [];
elseif prev_index < min(to_erase_indices) || prev_index > max(to_erase_indices)
    after_min_index = mod(min(to_erase_indices), numel(edge_arr))+1;
    before_max_index = mod(max(to_erase_indices)-2, numel(edge_arr))+1;
    while edge_arr(after_min_index).x == x && edge_arr(after_min_index).y == y
       after_min_index = mod(after_min_index, numel(edge_arr))+1;
    end
    while edge_arr(before_max_index).x == x && edge_arr(before_max_index).y == y
        before_max_index = mod(before_max_index-2, numel(edge_arr))+1;
    end
    
    option_a = edge_arr(after_min_index:before_max_index);
    option_b = [edge_arr(1:min(to_erase_indices)-1)...
        edge_arr(max(to_erase_indices)+1:end)];
    radius_a = max_radius_in_edge_list(option_a, N);
    radius_b = max_radius_in_edge_list(option_b, N);
    if radius_a > radius_b
        prev_index = before_max_index;
        prev_edge = edge_arr(prev_index);
        edge_arr=option_a;
        if numel(option_a)~= numel(option_b)
            grid = fill_points_with_3(grid, option_b);
        end
    else
        edge_arr=option_b;
        if numel(option_a)~= numel(option_b)
            grid = fill_points_with_3(grid, option_a);  
        end
    end
else
%     grid(prev_edge.x, prev_edge.y) = 3;
    next_index = mod(index, numel(edge_arr))+1;
    if max(to_erase_indices)< numel(edge_arr)-1
        next_index = max(to_erase_indices)+1;
    end
    next_edge = edge_arr(next_index);
    next_circled_to_start = false;
    while next_edge.x == x && next_edge.y == y
        new_next_index = mod(next_index, numel(edge_arr))+1;
        if new_next_index < next_index
            next_circled_to_start = true;
        end
        next_index = new_next_index;
        next_edge = edge_arr(next_index);
    end
    if next_circled_to_start || next_index == 2
        edge_arr = [edge_arr(next_index:end) edge_arr(1:next_index-1)]; % cycle the array so it's easy to work with
        to_erase_fn = @(i) edge_arr(i).x == x && edge_arr(i).y == y;
        to_erase_indices = find(arrayfun(to_erase_fn, 1:numel(edge_arr)));
        prev_index = prev_index - next_index+1;
        next_index=1;
    end
    if next_index < prev_index && next_index > min(to_erase_indices)
        edge_arr = edge_arr(next_index: prev_index);
%         grid(prev_edge.x, prev_edge.y) = 0;
    else
        before_min_index = mod(min(to_erase_indices)-2, numel(edge_arr))+1;
        while edge_arr(before_min_index).x == x && edge_arr(before_min_index).y == y
            before_min_index = mod(before_min_index-2, numel(edge_arr))+1;
        end
        after_min_index = mod(min(to_erase_indices), numel(edge_arr))+1;
        while edge_arr(after_min_index).x == x && edge_arr(after_min_index).y == y
            after_min_index = mod(after_min_index, numel(edge_arr))+1;
        end
        %         ax = edge_arr(after_min_index).x; ay= edge_arr(after_min_index).y;
%         contains3 =  @(i) edge_arr(i).x == prev_edge.x && edge_arr(i).y == prev_edge.y;
%         if edge_arr(after_min_index).x == prev_edge.x && edge_arr(after_min_index).y == prev_edge.y ...
%                 || edge_arr(before_min_index).x == edge_arr(next_index).x && ...
%                 edge_arr(before_min_index).y == edge_arr(next_index).y
%             %any(arrayfun(contains3, after_min_index : prev_index))
% %              grid(prev_edge.x, prev_edge.y) = 0;
%             %grid(ax,ay) == 3
%             edge_arr = [edge_arr(1:before_min_index), edge_arr(next_index:end)];
%             prev_edge = edge_arr(before_min_index);
%         else
%             edge_arr = edge_arr(after_min_index:prev_index);
            ax = edge_arr(after_min_index).x; ay = edge_arr(after_min_index).y;
            option_a = edge_arr(after_min_index:prev_index);
            if min(to_erase_indices) > 1 && next_index > 1
                option_b = [edge_arr(1:before_min_index), edge_arr(next_index:end)];
            else
                option_b = edge_arr(1:before_min_index);
            end
            radius_a = max_radius_in_edge_list(option_a, N);
            radius_b = max_radius_in_edge_list(option_b, N);
            if radius_a >= radius_b
                edge_arr = option_a;
                if numel(option_a)~= numel(option_b)
                    grid = fill_points_with_3(grid, option_b);
                end
            else
                edge_arr = option_b;
                prev_edge = edge_arr(before_min_index);
                if numel(option_a)~= numel(option_b)
                    grid = fill_points_with_3(grid, option_a);
                end
            end
%             grid(prev_edge.x, prev_edge.y) = 0;
%             grid(ax, ay) = 0;
%         end
    end
end
% if size_before - numel(edge_arr) > 12
%     stop = 'breakpoint'
% end
% edge_arr(tf) = [];

% then add new edges
new_index_of_prev_edge = @(i) edge_arr(i).x == prev_edge.x && ...
    edge_arr(i).y == prev_edge.y && edge_arr(i).from_x == prev_edge.from_x...
    && edge_arr(i).from_y == prev_edge.from_y ;
tf = arrayfun(new_index_of_prev_edge, 1:numel(edge_arr));
new_index = find(tf);
edge_arr = add_edges_after_index(edge_arr, new_index, grid, x, y, fr_x, fr_y);
end

function edge_arr = add_edges(edge_arr, index, grid)
% adds edges as if edge #index was chosen. does NOT remove any edge
y = edge_arr(index).y; fr_y = edge_arr(index).from_y;
x = edge_arr(index).x; fr_x = edge_arr(index).from_x;
p_x = x+1; m_x = x-1;
new_point_num = 1;
prev_index = mod(index - 2, numel(edge_arr)) +1;
if y > fr_y
    % if self.edge(prev_index).from_x < self.edge(min_index).from_x
    % go from '9' clocwise
    new_point =[];
    if grid(m_x, y) == 0
        new_point(new_point_num).x = m_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(m_x, y) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    if grid(x, y +1) == 0 %|| (grid(p_x, y) == 0 && grid(m_x, y) == 0)
        new_point(new_point_num).x = x;
        new_point(new_point_num).y = y+1;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(x, y +1) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    if grid(p_x, y) == 0
        new_point(new_point_num).x = p_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(p_x, y) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    edge_arr = [edge_arr(1:prev_index), new_point, ...
        edge_arr(prev_index + 1 : end)];
    
elseif fr_y == y %self.edge(prev_index).from_x == self.edge(min_index).from_x
    if  (fr_x-x == 1  && fr_x > x)% || (x == self.M && fr_x == 1) %self.edge(prev_index).from_y < self.edge(min_index).from_y
        % means we have something like
        %           :--|
        %           ...|
        %           x--|   <- min
        %           ...|   <- prev
        %           ___|
        % go from '6' clockwise
        new_point = [];
        if grid(x, y - 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y - 1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y - 1) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        if grid(m_x, y) == 0 %|| (grid(x, y - 1) == 0 && grid(x, y + 1) == 0)
            new_point(new_point_num).x = m_x;
            new_point(new_point_num).y = y;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(m_x, y) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        if grid(x, y + 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y+1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y + 1) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        edge_arr = [edge_arr(1:prev_index), new_point, ...
            edge_arr(prev_index + 1 : end)];
    else
        % means we have something like
        %           |--:
        %           |...    <- prev
        %           |--x    <- min
        %           |...
        %         __|___
        % go from '12' clockwise
        new_point_num = 1;
        new_point = [];
        if grid(x, y + 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y + 1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y + 1) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        if grid(p_x, y) == 0 %|| (grid(x, y + 1) == 0 && grid(x, y - 1) == 0)
            new_point(new_point_num).x = p_x;
            new_point(new_point_num).y = y;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(p_x, y) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        if grid(x, y - 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y - 1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y - 1) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        edge_arr = [edge_arr(1:prev_index), new_point, ...
            edge_arr(prev_index + 1 : end)];
    end
elseif y < fr_y
    % means we have something like
    %           |--|----
    %           |  x   ^
    %           |  ^   prev
    %           |..min
    %         __|___
    % go from '3' clockwise
    new_point_num = 1;
    new_point = [];
    if grid(p_x, y) == 0
        new_point(new_point_num).x = p_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(p_x, y) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    if grid(x, y - 1) == 0 %|| (grid(m_x, y) == 0 && grid(p_x, y) == 0)
        new_point(new_point_num).x = x;
        new_point(new_point_num).y = y - 1;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(x, y - 1) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    if grid(m_x, y) == 0
        new_point(new_point_num).x = m_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(m_x, y) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    edge_arr = [edge_arr(1:prev_index), new_point, ...
        edge_arr(prev_index + 1 : end)];
    
end
end

function edge_arr = add_edges_after_index(edge_arr, index, grid, x, y, ...
    fr_x, fr_y)
% adds edge after the edge with the given parameters (x, y, fr_x, fr_y) was
% chosen. adds the new edges at position #index.
% suitable for when edges were deleted before adding the new edges.
p_x = x+1; m_x = x-1;
new_point_num = 1;
% prev_index = mod(index - 2, numel(edge_arr)) +1;
if y > fr_y
    % if self.edge(prev_index).from_x < self.edge(min_index).from_x
    % go from '9' clocwise
    new_point =[];
    if grid(m_x, y) == 0
        new_point(new_point_num).x = m_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(m_x, y) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    if grid(x, y +1) == 0 %|| (grid(p_x, y) == 0 && grid(m_x, y) == 0)
        new_point(new_point_num).x = x;
        new_point(new_point_num).y = y+1;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(x, y +1) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    if grid(p_x, y) == 0
        new_point(new_point_num).x = p_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(p_x, y) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    edge_arr = [edge_arr(1:index), new_point, ...
        edge_arr(index + 1 : end)];
    
elseif fr_y == y %self.edge(prev_index).from_x == self.edge(min_index).from_x
    if  (fr_x-x == 1  && fr_x > x)% || (x == self.M && fr_x == 1) %self.edge(prev_index).from_y < self.edge(min_index).from_y
        % means we have something like
        %           :--|
        %           ...|
        %           x--|   <- min
        %           ...|   <- prev
        %           ___|
        % go from '6' clockwise
        new_point = [];
        if grid(x, y - 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y - 1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y - 1) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        if grid(m_x, y) == 0 %|| (grid(x, y - 1) == 0 && grid(x, y + 1) == 0)
            new_point(new_point_num).x = m_x;
            new_point(new_point_num).y = y;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(m_x, y) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        if grid(x, y + 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y+1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y + 1) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        edge_arr = [edge_arr(1:index), new_point, ...
            edge_arr(index + 1 : end)];
    else
        % means we have something like
        %           |--:
        %           |...    <- prev
        %           |--x    <- min
        %           |...
        %         __|___
        % go from '12' clockwise
        new_point_num = 1;
        new_point = [];
        if grid(x, y + 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y + 1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y + 1) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        if grid(p_x, y) == 0 %|| (grid(x, y + 1) == 0 && grid(x, y - 1) == 0)
            new_point(new_point_num).x = p_x;
            new_point(new_point_num).y = y;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(p_x, y) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        if grid(x, y - 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y - 1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y - 1) == 0;
            new_point(new_point_num).counter = 0;
            new_point_num = new_point_num + 1;
        end
        edge_arr = [edge_arr(1:index), new_point, ...
            edge_arr(index + 1 : end)];
    end
elseif y < fr_y
    % means we have something like
    %           |--|----
    %           |  x   ^
    %           |  ^   prev
    %           |..min
    %         __|___
    % go from '3' clockwise
    new_point_num = 1;
    new_point = [];
    if grid(p_x, y) == 0
        new_point(new_point_num).x = p_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(p_x, y) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    if grid(x, y - 1) == 0 %|| (grid(m_x, y) == 0 && grid(p_x, y) == 0)
        new_point(new_point_num).x = x;
        new_point(new_point_num).y = y - 1;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(x, y - 1) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    if grid(m_x, y) == 0
        new_point(new_point_num).x = m_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(m_x, y) == 0;
        new_point(new_point_num).counter = 0;
        new_point_num = new_point_num + 1;
    end
    edge_arr = [edge_arr(1:index), new_point, ...
        edge_arr(index + 1 : end)];
    
end
end

function edge = initial_edge(N)
x = N/2; y=N/2; p_x=N/2+1; m_x=N/2-1;
edge = struct('x', [], 'y', [], 'from_x', [], 'from_y', [], 'valid', []);
edge(1).x = m_x;
edge(1).y = y;
edge(1).from_x = x;
edge(1).from_y = y;
edge(1).valid = true;
edge(1).counter = 0;

edge(2).x = x;
edge(2).y = y+1;
edge(2).from_x = x;
edge(2).from_y = y;
edge(2).valid = true;
edge(2).counter = 0;

edge(3).x = p_x;
edge(3).y = y;
edge(3).from_x = x;
edge(3).from_y = y;
edge(3).valid = true;
edge(3).counter = 0;

edge(4).x = x;
edge(4).y = y-1;
edge(4).from_x = x;
edge(4).from_y = y;
edge(4).valid = true;
edge(4).counter = 0;
end

function [xnew,ynew] = make_step_vec(x, y, dones)
step = ceil(rand(length(x),1)*4); % choose direction for the walker
step(dones) = 0; % walkers that already hit the aggregate
xnew=x; ynew=y;
xnew(step==1)=x(step==1)-1;
xnew(step==2)=x(step==2)+1;
ynew(step==3)= ynew(step==3)-1;
ynew(step==4)= ynew(step==4)+1;
end

function [x, y] = check_radius_and_replace(x, y, rmax, dones, N, grid)
r_sq = x.*x + y.*y;
out = (r_sq >= 1.5*(rmax+2)^2);
% generate new walkers intead of the escaped ones
while any(out) % Loop until we get a good starting point for the walker
    theta = rand(length(x), 1)*2*pi;% Pick a random angle
    x(out) = ceil((rmax+2)*cos(theta(out))); % Starting x and y position
    y(out) = ceil((rmax+2)*sin(theta(out))); % of the walker
    if (all(abs(x(out)) < N/2) && all(abs(y(out)) < N/2))
        inds = sub2ind([N N], x(out) + N/2,y(out)+ N/2);
        if all(grid(inds)==0); break ; end
    end
end
end

function radius = max_radius_in_edge_list(edge_arr, N)
r_arr = arrayfun(@(i) (edge_arr(i).x-N/2)^2 + (edge_arr(i).y-N/2)^2,...
    1:numel(edge_arr));
r_arr = sqrt(r_arr);
radius = max(r_arr);
end

function grid=fill_points_with_3(grid, edge_arr)
for i=1:numel(edge_arr)
    x = edge_arr(i).x;
    y= edge_arr(i).y;
    grid(x,y)=3;
end
end


