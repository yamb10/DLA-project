function DLA_edge_walk(N, edge_walk_end_prob)

cc = [1 1 1; 0 0 1]; % Blue and yellow (RGB values)
colormap(cc); % Put these two colors into the color map
set(0,'defaultaxesfontsize',20) % Make the default font size bigger
grid = zeros(N,N); % Domain where growth can occur.
% Empty sites are 0, cluster sites 2, placeholders are 3
grid(N/2,N/2) = 2; % Seed of the cluster
edge = initial_edge(N);
rmax = 1; % Outer edge of the cluster
h = image(grid);
axis image
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
% Keep the cluster growing until it touches any edge of the
% simulation cell.
while (all(grid(:,1)==0) && all(grid(:,N)==0) && ...
        all(grid(1,:)==0) && all(grid(N,:)==0))
    r = rmax+10; % Pick a distance just beyond the edge of the cluster
    if r>N/2; break; end % Stop when the cluster gets too big
    while 1==1 % Loop until we get a good starting point for the walker
        theta = rand()*2*pi;% Pick a random angle
        x = ceil(r*cos(theta)); % Starting x and y position
        y = ceil(r*sin(theta)); % of the walker
        if (abs(x) < N/2 && abs(y) < N/2)
            if grid(x+N/2,y+N/2)==0; break ; end
        end
    end
    % Move the walker randomly left, right, up, or down
    % until it 'escapes' or is stuck to the cluster
    while 1==1
        step = ceil(rand()*4); % choose direction for the walker
        if step == 1
            xnew = x - 1; ynew = y;
        elseif step == 2
            xnew = x + 1; ynew = y;
        elseif step == 3
            ynew = y - 1; xnew = x;
        else
            ynew = y + 1; xnew = x;
        end
        if xnew^2+ynew^2 > (rmax+10)^2*1.5; break; end % The walker 'escaped'
        % If the walker is outside of the box but not too far away, keep going
%         if (abs(xnew)>=N/2 || abs(ynew)>=N/2); continue; end
        % Otherwise, check to see if it hit the cluster
        if ~(abs(xnew)>=N/2 || abs(ynew)>=N/2) && ...
                grid(xnew+N/2,ynew+N/2)==2 % Then the walker 'stuck' to the cluster
            [grid, edge] = insert_point(grid, edge, x+N/2, y+N/2,...
                xnew+N/2, ynew+N/2, edge_walk_end_prob); % Stick the walker at its previous step
            if (x^2+y^2)>rmax^2; rmax = sqrt(x^2+y^2); end
            break % Break out of the loop and start a new walker
        end
        x = xnew; y = ynew; % Update the walker position if it hasn't stuck yet
    end
    set(h,'cdata',grid)
    drawnow
end
end


function [grid, edge_arr] = insert_point(grid, edge_arr, x, y, fr_x, fr_y, prob)
    new_edge = choose_new_edge(edge_arr, x, y, fr_x, fr_y, prob);
    
    new_x = new_edge.x; new_y = new_edge.y;
    if grid(new_x, new_y) == 2
        error = 'dk'
    end
    grid(new_x, new_y) = 2;
    
    edge_arr = update_edge_and_erase_holes(edge_arr, new_edge, grid);
end

function [new_edge] = choose_new_edge(edge_arr, x, y, fr_x, fr_y, prob)
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
        if step_size > 0
            step_size = step_size -1;
        else
            step_size = step_size +1;
        end
        new_edge_index = mod(index - 1 + step_size, numel(edge_arr)) + 1;
        new_edge = edge_arr(new_edge_index);
    end
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

function [edge_arr] = update_edge_and_erase_holes(edge_arr, new_edge, grid)
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
if prev_index < min(to_erase_indices) || prev_index > max(to_erase_indices)
    edge_arr(min(to_erase_indices):max(to_erase_indices)) = [];
else
    grid(prev_edge.x, prev_edge.y) = 3;
    next_index = mod(index, numel(edge_arr))+1;
    next_edge = edge_arr(next_index);
    while next_edge.x == x && next_edge.y == y
        next_index = mod(next_index, numel(edge_arr))+1;
        next_edge = edge_arr(next_index);
    end
    if next_index < prev_index
        edge_arr = edge_arr(next_index: prev_index);
    else
        before_min_index = mod(min(to_erase_indices)-2, numel(edge_arr))+1;
        after_min_index = mod(min(to_erase_indices), numel(edge_arr))+1;
        while edge_arr(after_min_index).x == x && edge_arr(after_min_index).y == y
            after_min_index = mod(after_min_index, numel(edge_arr))+1;
        end
%         ax = edge_arr(after_min_index).x; ay= edge_arr(after_min_index).y;
        contains3 =  @(i) edge_arr(i).x == prev_edge.x && edge_arr(i).y == prev_edge.y;
        if any(arrayfun(contains3, after_min_index : prev_index))
            %grid(ax,ay) == 3
            edge_arr = [edge_arr(1:before_min_index), edge_arr(next_index:end)];
            prev_edge = edge_arr(before_min_index);
        else
            edge_arr = [edge_arr(after_min_index:prev_index)];
        end
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
        new_point_num = new_point_num + 1;
    end
    if grid(x, y +1) == 0 %|| (grid(p_x, y) == 0 && grid(m_x, y) == 0)
        new_point(new_point_num).x = x;
        new_point(new_point_num).y = y+1;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(x, y +1) == 0;
        new_point_num = new_point_num + 1;
    end
    if grid(p_x, y) == 0
        new_point(new_point_num).x = p_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(p_x, y) == 0;
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
            new_point_num = new_point_num + 1;
        end
        if grid(m_x, y) == 0 %|| (grid(x, y - 1) == 0 && grid(x, y + 1) == 0)
            new_point(new_point_num).x = m_x;
            new_point(new_point_num).y = y;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(m_x, y) == 0;
            new_point_num = new_point_num + 1;
        end
        if grid(x, y + 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y+1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y + 1) == 0;
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
            new_point_num = new_point_num + 1;
        end
        if grid(p_x, y) == 0 %|| (grid(x, y + 1) == 0 && grid(x, y - 1) == 0)
            new_point(new_point_num).x = p_x;
            new_point(new_point_num).y = y;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(p_x, y) == 0;
            new_point_num = new_point_num + 1;
        end
        if grid(x, y - 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y - 1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y - 1) == 0;
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
        new_point_num = new_point_num + 1;
    end
    if grid(x, y - 1) == 0 %|| (grid(m_x, y) == 0 && grid(p_x, y) == 0)
        new_point(new_point_num).x = x;
        new_point(new_point_num).y = y - 1;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(x, y - 1) == 0;
        new_point_num = new_point_num + 1;
    end
    if grid(m_x, y) == 0
        new_point(new_point_num).x = m_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(m_x, y) == 0;
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
        new_point_num = new_point_num + 1;
    end
    if grid(x, y +1) == 0 %|| (grid(p_x, y) == 0 && grid(m_x, y) == 0)
        new_point(new_point_num).x = x;
        new_point(new_point_num).y = y+1;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(x, y +1) == 0;
        new_point_num = new_point_num + 1;
    end
    if grid(p_x, y) == 0
        new_point(new_point_num).x = p_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(p_x, y) == 0;
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
            new_point_num = new_point_num + 1;
        end
        if grid(m_x, y) == 0 %|| (grid(x, y - 1) == 0 && grid(x, y + 1) == 0)
            new_point(new_point_num).x = m_x;
            new_point(new_point_num).y = y;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(m_x, y) == 0;
            new_point_num = new_point_num + 1;
        end
        if grid(x, y + 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y+1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y + 1) == 0;
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
            new_point_num = new_point_num + 1;
        end
        if grid(p_x, y) == 0 %|| (grid(x, y + 1) == 0 && grid(x, y - 1) == 0)
            new_point(new_point_num).x = p_x;
            new_point(new_point_num).y = y;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(p_x, y) == 0;
            new_point_num = new_point_num + 1;
        end
        if grid(x, y - 1) == 0
            new_point(new_point_num).x = x;
            new_point(new_point_num).y = y - 1;
            new_point(new_point_num).from_x = x;
            new_point(new_point_num).from_y = y;
            new_point(new_point_num).valid = grid(x, y - 1) == 0;
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
        new_point_num = new_point_num + 1;
    end
    if grid(x, y - 1) == 0 %|| (grid(m_x, y) == 0 && grid(p_x, y) == 0)
        new_point(new_point_num).x = x;
        new_point(new_point_num).y = y - 1;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(x, y - 1) == 0;
        new_point_num = new_point_num + 1;
    end
    if grid(m_x, y) == 0
        new_point(new_point_num).x = m_x;
        new_point(new_point_num).y = y;
        new_point(new_point_num).from_x = x;
        new_point(new_point_num).from_y = y;
        new_point(new_point_num).valid = grid(m_x, y) == 0;
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

edge(2).x = x;
edge(2).y = y+1;
edge(2).from_x = x;
edge(2).from_y = y;
edge(2).valid = true;

edge(3).x = p_x;
edge(3).y = y;
edge(3).from_x = x;
edge(3).from_y = y;
edge(3).valid = true;

edge(4).x = x;
edge(4).y = y-1;
edge(4).from_x = x;
edge(4).from_y = y;
edge(4).valid = true;

end


