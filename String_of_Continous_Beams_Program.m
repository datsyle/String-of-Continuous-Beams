%% 
% Inputting and organizing information to run the String of Continous Beams 
% code:


% This is the format that the inputdata should be placed in to properly run
% the String of Continous Beams Code


% the way information in this matrix is arranged is that the rows above and below
% each other in the matrix  are paired together to represent the positions in [x y z] at the beginning and ends 
% of each segment in the system
 


system_filament_segment_positioning_data = [  0 10 0 ;
                                              0 9 0  ;
                                              0 9 0  ;
                                              1 8 3  ;
                                              0 9 0  ;
                                             -1 8 4  ;
                                              0 9 0  ;
                                              0 8 -3 ]

[total_number_of_segment_endings, ~ ] = size(system_filament_segment_positioning_data)



% creating separate matrices for the positions of the front and end nodes in [x y z] for each
% segment input variable "system_filament_segment_positioning_data"



front_node_positions_of_segments = [];
end_node_positions_of_segments = [];

for row_number = 1:total_number_of_segment_endings

    is_odd = rem(row_number, 2) == 1;

    if is_odd == true
        front_node_positions_of_segments = [front_node_positions_of_segments ; system_filament_segment_positioning_data(row_number,:)];
    
    else 
        end_node_positions_of_segments = [end_node_positions_of_segments ; system_filament_segment_positioning_data(row_number,:)];

    end

end

front_node_positions_of_segments

end_node_positions_of_segments



% calculating the total amount of nodes in the whole system



total_number_of_nodes = 0;

for number_of_front_or_end_nodes = 1: (total_number_of_segment_endings / 2)

    if front_node_positions_of_segments(number_of_front_or_end_nodes, :) == end_node_positions_of_segments(number_of_front_or_end_nodes, :)
        total_number_of_nodes = total_number_of_nodes;

    else
        total_number_of_nodes = total_number_of_nodes + 1;

    end

end

total_number_of_nodes



% separating the original input data into the seperate position vectors x y z
% and plotting them out to see the initial system

x = system_filament_segment_positioning_data(: , 1);

y = system_filament_segment_positioning_data(: , 2);

z = system_filament_segment_positioning_data(: , 3);

x_start = front_node_positions_of_segments(: , 1);

y_start = front_node_positions_of_segments(: , 2);

z_start = front_node_positions_of_segments(: , 3);

x_end = end_node_positions_of_segments(: , 1);

y_end = end_node_positions_of_segments(: , 2);

z_end = end_node_positions_of_segments(: , 3);

plot3(x, z, y)

% note : the reason why the y and z are swapped in the scatterplot function is
    %   because MatLab's scatter graph assumes z to be in the upwards direction,
    %   while y is the assumed upwards direction for the String of Continous
    %   Beams Code
%% 
% Solving the SOCB Code:
% 
% Step 1: Solving for reaction forces at every node

% it is assumed that the forces and moments applied to the frame created
% will be at the first "front" node of the first segment created in the
% code


% Node at which the forces and moments will be applied
NP = system_filament_segment_positioning_data(1 , :);

% Forces in the x-direction to be applied
Fx = 0

% Forces in the y-direction to be applied
Fy = -10

% Forces in the z-direction to be applied
Fz = 0

% Moments in the xy-plane to be applied
Mxy = 0

% Moments in the yz-plane to be applied
Myz = 0

% Moments in the xz-plane to be applied (twisting)
Mxz = 0


%===============================================================================%
% This is to check if NP is a valid position to put forces and moments


NP_is_valid = true;
for checking_if_NP_is_valid = 1:total_number_of_segment_endings

    if NP == system_filament_segment_positioning_data(checking_if_NP_is_valid , :);

        NP = true;
        break

    else

        NP_is_valid = false;

    end

end

if NP_is_valid == false
    invalid = MException('MyComponent:noSuchVariable','Position chosen for variable NP is invalid %s ');
    throw(invalid)
end



%===============================================================================%

% calculating the lengths of each segment

length_of_segments = [];

for segment_number = 1:total_number_of_nodes
    length = sqrt(sum((front_node_positions_of_segments(segment_number, :) - end_node_positions_of_segments(segment_number, :)).^2));
    length_of_segments = [length_of_segments ; length];
end

length_of_segments

length_of_segments_x = front_node_positions_of_segments(:, 1) - end_node_positions_of_segments(: , 1)

length_of_segments_y = front_node_positions_of_segments(:, 2) - end_node_positions_of_segments(: , 2)

length_of_segments_z = front_node_positions_of_segments(:, 3) - end_node_positions_of_segments(: , 3)



% Solving for Reaction Forces in the frame:


% This is done by solving for which members are connected, turning them into unit vectors and creating
% turning them into a linear system which will solve for the reaction forces in
% terms of moments and forces


% turning the members into equivalent unit vectors



unit_vector_member_data = [];

for number_of_node = 1:total_number_of_nodes

    row_unit_vector_members = [  end_node_positions_of_segments(number_of_node , :) - front_node_positions_of_segments(number_of_node , :) ] ./length_of_segments(number_of_node , :);

    unit_vector_member_data = [unit_vector_member_data ; row_unit_vector_members];

end

unit_vector_member_data



% Finding how many reaction forces there will be and assigning values to be
% input to stiffness matrix

total_number_of_connections = [];

for end_nodes = 1:total_number_of_nodes

    number_of_connections = 0;

    for front_nodes = 1:total_number_of_nodes

        if end_node_positions_of_segments(end_nodes , :) == front_node_positions_of_segments(front_nodes , :)

            number_of_connections = number_of_connections + 1;

        end

    end

    if number_of_connections == 0

        number_of_connections = [];
    end

    total_number_of_connections = [total_number_of_connections ; number_of_connections];

end

total_number_of_connections



% This program will support a maximum of 3 connected members connected to
% the end of another member


% this is what the following stiffness matrix looks like on paper, where the maxumium number of members, 3, are attached to the end of another member:

% member_1_unit_vector_x   member_2_unit_vector_x   member_3_unit_vector_x          0                       0                        0            
% member_1_unit_vector_y   member_2_unit_vector_y   member_3_unit_vector_y          0                       0                        0 
% member_1_unit_vector_z   member_2_unit_vector_z   member_3_unit_vector_z          0                       0                        0
%           0                         0                       0            member_1_unit_vector_z   member_2_unit_vector_z   member_3_unit_vector_z
%           0                         0                       0            member_1_unit_vector_x   member_2_unit_vector_x   member_3_unit_vector_x
%           0                         0                       0            member_1_unit_vector_y   member_2_unit_vector_y   member_3_unit_vector_y


% This is for the first beam that is solved in the program, the variables
% used in repetition will be different from the first beam solved.

forces_x = -1 * Fx

forces_y = -1 *  Fy

forces_z = -1 * Fz

moments_xy = -1 * Mxy

moments_yz = -1 * Myz

moments_xz = -1 * Mxz

[height_total_forces_and_moments_xyz,~] = size(end_node_positions_of_segments)

total_forces_and_moments_xyz = zeros(height_total_forces_and_moments_xyz, 6) 

[number_of_iterations, ~] = size(total_number_of_connections)

% The next_line variable is to decide how many lines to move down and take data from when iterating
% through unit_vector_member_data array based on the amount of connections
% a node has and adds it to n, which will determine where the code will
% start to read the information from, so it doesn't double count

n = 1

for i = 1:number_of_iterations

    if total_number_of_connections(i, 1) == 1

        next_line = 1
    
    elseif total_number_of_connections(i, 1) == 2

        next_line = 2
    
    elseif total_number_of_connections(i, 1) == 3

        next_line = 3
    
    end

    applied_forces_input = [forces_x ;
                            forces_y ;
                            forces_z ;
                            moments_xy ;
                            moments_yz ;
                            moments_xz ] ;

    % This forces and moment adjustment below is to satisfy equilibrium
    % between applied forces created from applied moments and vice versa


    applied_forces_adjustment = [- (moments_xz ./ length_of_segments_z(n , 1)) ;
                                 - (moments_xy ./ length_of_segments_x(n , 1)) ;
                                 - (moments_yz ./ length_of_segments_y(n , 1)) ;
                                 + (forces_y * length_of_segments_z(n , 1) - forces_z * length_of_segments_y(n , 1)) ;
                                 + (forces_z * length_of_segments_x(n , 1) - forces_x * length_of_segments_z(n , 1)) ;
                                 + (forces_x * length_of_segments_y(n , 1) - forces_y * length_of_segments_x(n , 1))] ;

    applied_forces_adjustment(isnan(applied_forces_adjustment)) = 0 ;

    applied_forces = applied_forces_input + applied_forces_adjustment

    stiffness_matrix = zeros(6, 6);




    if total_number_of_connections(i , 1) == 1

        stiffness_matrix(1, 1) = unit_vector_member_data(n + 1, 1);

        stiffness_matrix(2, 1) = unit_vector_member_data(n + 1, 2);

        stiffness_matrix(3, 1) = unit_vector_member_data(n + 1, 3);

        stiffness_matrix(4, 4) = unit_vector_member_data(n + 1, 3);

        stiffness_matrix(5, 4) = unit_vector_member_data(n + 1, 1);

        stiffness_matrix(6, 4) = unit_vector_member_data(n + 1, 2)

    elseif total_number_of_connections(i, 1) == 2

        stiffness_matrix(1:3, 1) = unit_vector_member_data(n + 1 ,:);

        stiffness_matrix(1:3, 2) = unit_vector_member_data(n + 2, :);

        stiffness_matrix(4, 4:5) = unit_vector_member_data(n + 1: n + 2, 3);

        stiffness_matrix(5, 4:5) = unit_vector_member_data(n + 1: n + 2, 1)

    elseif total_number_of_connections(i, 1) == 3

        stiffness_matrix(1:3, 1) = unit_vector_member_data(n + 1 ,:);

        stiffness_matrix(1:3, 2) = unit_vector_member_data(n + 2, :);

        stiffness_matrix(1:3, 3) = unit_vector_member_data(n + 3, :);

        stiffness_matrix(4, 4:6) = unit_vector_member_data(n + 1: n + 3, 1);

        stiffness_matrix(5, 4:6) = unit_vector_member_data(n + 1: n + 3, 2);

        stiffness_matrix(6, 4:6) = unit_vector_member_data(n + 1: n + 3, 3);

        stiffness_matrix(4, 1:3) = unit_vector_member_data(n + 1: n + 3, 2) .* length_of_segments_z(n + 1, 1) - unit_vector_member_data(n + 1: n + 3, 3) .* length_of_segments_y(n + 1, 1);

        stiffness_matrix(5, 1:3) = unit_vector_member_data(n + 1: n + 3, 3) .* length_of_segments_x(n + 1, 1) - unit_vector_member_data(n + 1: n + 3, 1) .* length_of_segments_z(n + 1, 1);
        
        stiffness_matrix(6, 1:3) = unit_vector_member_data(n + 1: n + 3, 1) .* length_of_segments_y(n + 1, 1) - unit_vector_member_data(n + 1: n + 3, 2) .* length_of_segments_x(n + 1, 1)


        stiffness_matrix(isnan(stiffness_matrix)) = 0

        

    end



    n = n + next_line

    forces_and_moment_reactions = stiffness_matrix \ applied_forces 

    sum_forces_x =   - forces_and_moment_reactions(1, 1)

    sum_forces_y =   - forces_and_moment_reactions(2, 1)

    sum_forces_z =   - forces_and_moment_reactions(3, 1)

    sum_moments_member1 = - forces_and_moment_reactions(4, 1)

    sum_moments_member2 = - forces_and_moment_reactions(5, 1)

    sum_moments_member3 = - forces_and_moment_reactions(6, 1)

    nodal_reactions = [sum_forces_x, sum_forces_y, sum_forces_z, sum_moments_member1, sum_moments_member2, sum_moments_member3];



    if total_number_of_connections(i , 1) == 1

        total_forces_and_moments_xyz(  i  , 1:6) = -1* applied_forces

        total_forces_and_moments_xyz(i + 1, 1:3) = nodal_reactions(1, 1) .* unit_vector_member_data(i + 1, :);

        total_forces_and_moments_xyz(i + 1, 4:6) = nodal_reactions(1, 4) .* unit_vector_member_data(i + 1, :)

    elseif total_number_of_connections(i, 1) == 2

        total_forces_and_moments_xyz(  i  , 1:6) = -1* applied_forces

        total_forces_and_moments_xyz(i + 1, 1:3) = nodal_reactions(1, 1) .* unit_vector_member_data(i + 1, :);

        total_forces_and_moments_xyz(i + 2, 1:3) = nodal_reactions(1, 2) .* unit_vector_member_data(i + 2, :);

        total_forces_and_moments_xyz(i + 1, 4:6) = nodal_reactions(1, 4) .* unit_vector_member_data(i + 1, :);

        total_forces_and_moments_xyz(i + 2, 4:6) = nodal_reactions(1, 5) .* unit_vector_member_data(i + 2, :)
        
    elseif total_number_of_connections(i, 1) == 3

        total_forces_and_moments_xyz(  i  , 1:6) = -1* applied_forces

        total_forces_and_moments_xyz(i + 1, 1:3) = nodal_reactions(1, 1) .* unit_vector_member_data(i + 1, :);

        total_forces_and_moments_xyz(i + 2, 1:3) = nodal_reactions(1, 2) .* unit_vector_member_data(i + 2, :);

        total_forces_and_moments_xyz(i + 3, 1:3) = nodal_reactions(1, 3) .* unit_vector_member_data(i + 3, :);

        total_forces_and_moments_xyz(i + 1, 4:6) = nodal_reactions(1, 4) .* unit_vector_member_data(i + 1, :);

        total_forces_and_moments_xyz(i + 2, 4:6) = nodal_reactions(1, 5) .* unit_vector_member_data(i + 2, :);

        total_forces_and_moments_xyz(i + 3, 4:6) = nodal_reactions(1, 6) .* unit_vector_member_data(i + 3, :)

    end

    

end

% total_forces_and_moments_xyz is organized as follows:
%                | Forces in x | Forces in y | Forces in z | Moment around x-axis (or on yz - plane) | Moment around y-axis (or on xz - plane)| Moment around z-axis (or on xy-plane) |
%   Member 1 :   |             |             |             |                                         |                                        |                                       |  
%   Member 2 :   |             |             |             |                                         |                                        |                                       |
%   Member ... : |             |             |             |                                         |                                        |                                       |


% note to self: this solver solves the reaction forces in the xyz
% coordinate system and then should translate itself into the spherical
% coordinates system

% translating xyz forces and moments into spherical coordinate forces and
% moments

rho = [];
phi = [];
theta = [];

for number_of_member = 1: total_number_of_nodes   
    new_rho = sqrt( (x_start(number_of_member, :) - x_end(number_of_member, :)).^2 + (y_start(number_of_member, :) - y_end(number_of_member, :)).^2 + (z_start(number_of_member, :) - z_end(number_of_member, :)).^2 );
    rho = [rho; -new_rho];
end

rho

for number_of_member = 1: total_number_of_nodes   
    new_theta = acos((z_start(number_of_member, :) - z_end(number_of_member, :)) ./ rho(number_of_member, :));
    theta = [theta; -new_theta];
end

theta

for number_of_member = 1: total_number_of_nodes   
    new_phi = asin( (y_start(number_of_member, :) - y_end(number_of_member, :)) ./ ( rho(number_of_member, :) .* sin(theta(number_of_member, :)) ) );
    phi = [phi; -new_phi];
end

phi

total_forces_and_moments_spherical = zeros(height_total_forces_and_moments_xyz, 6);

for i = 1:number_of_iterations


    xyz_to_spherical_transformation_matrix = zeros(3,3);

    xyz_to_spherical_transformation_matrix(1, 1) = sin(theta(i + 1, 1))*cos(phi(i + 1, 1));

    xyz_to_spherical_transformation_matrix(1, 2) = sin(theta(i + 1, 1))*sin(phi(i + 1, 1));

    xyz_to_spherical_transformation_matrix(1, 3) = cos(theta(i + 1, 1));

    xyz_to_spherical_transformation_matrix(2, 1) = cos(theta(i + 1, 1))*cos(phi(i + 1, 1));

    xyz_to_spherical_transformation_matrix(2, 2) = cos(theta(i + 1, 1))*sin(phi(i + 1, 1));

    xyz_to_spherical_transformation_matrix(2, 3) = -1*sin(theta(i + 1, 1));

    xyz_to_spherical_transformation_matrix(3, 1) = -1*sin(phi(i + 1, 1));

    xyz_to_spherical_transformation_matrix(3, 2) = cos(phi(i + 1, 1));

    xyz_to_spherical_transformation_matrix(3, 3) = 1;


    if total_number_of_connections(i , 1) == 1

        new_spherical_forces = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(  i  , 1:3).');

        total_forces_and_moments_spherical(  i  , 1:3) =  new_spherical_forces.';

        new_spherical_forces = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 1, 1:3).');

        total_forces_and_moments_spherical(i + 1, 1:3) =  new_spherical_forces.';

        new_spherical_moments = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(  i  , 4:6).');

        total_forces_and_moments_spherical(  i  , 4:6) =  new_spherical_moments.';

        new_spherical_moments = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 1, 4:6).');

        total_forces_and_moments_spherical(i + 1, 4:6) =  new_spherical_moments.';

    elseif total_number_of_connections(i, 1) == 2

        new_spherical_forces = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(  i  , 1:3).');

        total_forces_and_moments_spherical(  i  , 1:3) =  new_spherical_forces.';
        
        new_spherical_forces = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 1, 1:3).');

        total_forces_and_moments_spherical(i + 1, 1:3) =  new_spherical_forces.';

        new_spherical_forces = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 2, 1:3).');

        total_forces_and_moments_spherical(i + 2, 1:3) =  new_spherical_forces.';

        new_spherical_moments = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(  i  , 4:6).');

        total_forces_and_moments_spherical(  i  , 4:6) =  new_spherical_moments.';

        new_spherical_moments = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 1, 4:6).');

        total_forces_and_moments_spherical(i + 1, 4:6) =  new_spherical_moments.';

        new_spherical_moments = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 2, 4:6).');

        total_forces_and_moments_spherical(i + 2, 4:6) =  new_spherical_moments.';


    elseif total_number_of_connections(i, 1) == 3

        new_spherical_forces = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(  i  , 1:3).');

        total_forces_and_moments_spherical(  i  , 1:3) =  new_spherical_forces.';

        new_spherical_forces = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 1, 1:3).');

        total_forces_and_moments_spherical(i + 1, 1:3) =  new_spherical_forces.';

        new_spherical_forces = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 2, 1:3).');

        total_forces_and_moments_spherical(i + 2, 1:3) =  new_spherical_forces.';

        new_spherical_forces = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 3, 1:3).');

        total_forces_and_moments_spherical(i + 3, 1:3) =  new_spherical_forces.';

        new_spherical_moments = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(  i  , 4:6).');

        total_forces_and_moments_spherical(  i  , 4:6) =  new_spherical_moments.';

        new_spherical_moments = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 1, 4:6).');

        total_forces_and_moments_spherical(i + 1, 4:6) =  new_spherical_moments.';

        new_spherical_moments = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 2, 4:6).');

        total_forces_and_moments_spherical(i + 2, 4:6) =  new_spherical_moments.';

        new_spherical_moments = xyz_to_spherical_transformation_matrix * (total_forces_and_moments_xyz(i + 3, 4:6).');

        total_forces_and_moments_spherical(i + 3, 4:6) =  new_spherical_moments.';

      
    end

end

total_forces_and_moments_spherical
    

% total_forces_and_moments_spherical is organized as follows:
%                | Forces in r | Forces in theta | Forces in phi | Moment around r-axis (or around r-axis) | Moment in theta  (idk what this means)| Moment around in phi |
%   Member 1 :   |             |                 |               |                                         |                                       |                      |  
%   Member 2 :   |             |                 |               |                                         |                                       |                      |
%   Member ... : |             |                 |               |                                         |                                       |                      |




% applying the coupling matrix to the system:

% length_span variable is what will be used to calculate the integrations
% as it will be the variable length for a given segment




a_theta = 0;

b_theta = 0;

c_theta = 0;

d_theta = 0;

e_theta = 0;

f_theta = 0;

a_phi = 0;

b_phi = 0;

c_phi = 0;

d_phi = 0;

e_phi = 0;

f_phi = 0;


height_of_contour_equations_array = height_total_forces_and_moments_xyz * 2;

contour_equation_coefficients_and_steps = zeros(height_of_contour_equations_array, 7);

for i = 1:total_number_of_nodes

    step_size = length_of_segments(i , 1) / ( 10 );

    contour_equation_coefficients_and_steps(i : i + 1 , 7 ) = step_size;

    length_span = 0 : step_size : length_of_segments(i , 1);

    unknown_coefficients_theta = [ a_theta , 
                                   b_theta ,
                                   c_theta ,
                                   d_theta ,
                                   e_theta ,
                                   f_theta ];

    unknown_coefficients_theta = [ a_phi , 
                                   b_phi ,
                                   c_phi ,
                                   d_phi ,
                                   e_phi ,
                                   f_phi ];

    coupling_matrix = [ 1                0                            0                                      0                       0                              0                  ;
                        0                1                            0                                      0                       0                              0                  ;
                        0                0                            1                                      0                       0                              0                  ;
                        1       length_of_segments(i , 1)     length_of_segments(i, 1)^ 2     length_of_segments(i, 1)^3     length_of_segments(i, 1)^4     length_of_segments(i, 1)^5 ;
                        0                1                    2*length_of_segments(i, 1)    3*length_of_segments(i, 1)^2   4*length_of_segments(i, 1)^3   5*length_of_segments(i, 1)^4 ;
                        0                0                            2                     6*length_of_segments(i, 1)    12*length_of_segments(i, 1)^2  20*length_of_segments(i, 1)^3 ];

    input_forces_and_moments = [                 theta(i, 1)                ,
                                  total_forces_and_moments_spherical(i , 2) ,
                                  total_forces_and_moments_spherical(i , 5) ,
                                             theta(i, 1)                    ,
                                 -total_forces_and_moments_spherical(i , 2) ,
                                 -total_forces_and_moments_spherical(i , 5) ];

    unknown_coefficients_theta = input_forces_and_moments \ coupling_matrix; 

    input_forces_and_moments = [                 phi(i, 1)                  ,
                                  total_forces_and_moments_spherical(i , 3) ,
                                  total_forces_and_moments_spherical(i , 6) ,
                                             phi(i, 1)                      ,
                                 -total_forces_and_moments_spherical(i , 3) ,
                                 -total_forces_and_moments_spherical(i , 6) ];

    unknown_coefficients_phi = input_forces_and_moments \ coupling_matrix; 

    contour_equation_coefficients_and_steps(i + (i - 1) , 1:6) = unknown_coefficients_theta;

    contour_equation_coefficients_and_steps(    2 * i   , 1:6) = unknown_coefficients_phi;


end    

contour_equation_coefficients_and_steps


% the information in variable contour_equation_coefficent_and_steps is
% organized as follows:

% |                                  | a_coefficient | b_coefficient | c_coefficient | d_coefficient | e_coefficient |f_coefficient|step intervals over length|
% | Member 1 theta contour equation: |               |               |               |               |               |             |                          |
% | Member 1 phi contour equation:   |               |               |               |               |               |             |                          |
% | Member 2 theta contour equation: |               |               |               |               |               |             |                          |
% | Member 2 phi contour equation:   |               |               |               |               |               |             |                          |

% Using these coefficients, the displacements can then be plotted using the
% following equation:

% theta(s) = a_coefficient_theta + b_coefficient_theta * s + c_coefficient_theta * s^2 + d_coefficient_theta * s^3 + e_coefficient_theta * s^4 + f_coefficient_theta * s^5

% phi(s) = a_coefficient_phi + b_coefficient_phi * s + c_coefficient_phi * s^2 + d_coefficient_phi * s^3 + e_coefficient_phi * s^4 + f_coefficient_phi * s^5