%CONSTS
MIN_RADIUS = 10;
DELTA_R = 2;
PARTICLE_SYMB = 2;


function [x, y, fr_x, fr_y] = simulate_single_particle(grid, ...
    starting_radius)
    N = length(grid);
    x, y = random_starting_point(starting_radius,N);
    new_x, new_y = make_step(x, y);
    while (~x_y_in_grid(new_x, new_y)) || grid(new_x, new_y) ~= PARTICLE_SYMB
        x = new_x;
        y = new_y;
        new_x, new_y = make_step(x, y);
    end
    fr_x = new_x; fr_y = new_y;
end

function [x, y] = random_starting_point(r, N)
    theta = pi * rand();
    R = r;
    x = floor(cos(theta) * R) + N/2;
    y = floor(sin(theta) * R) + N/2;
end

function [new_x, new_y] = make_step(x, y)
    direction = randi(4);
    switch direction
        case 1
            new_x = x + 1;
            new_y = y;
        case 2
            new_x = x - 1;
            new_y = y;
        case 3
            new_x = x;
            new_y = y + 1;
        case 4
            new_x = x;
            new_y = y - 1;
    end
end

function res = x_y_in_grid(x, y, grid)
    N = length(grid);
    res = (x > 0) && (x <= N) && (y > 0) && (y <= N);
end