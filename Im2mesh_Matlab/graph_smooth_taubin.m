function Vout = graph_smooth_taubin(Vin, F, feat, ntype, geom, lamda, mu, num_iters)
% graph_smooth_taubin: Smooth selected vertices with constrained Taubin
% smoothing.
%
% This function treats the selected part of a triangular surface mesh as a
% graph. Each movable vertex is connected to neighboring vertices by graph
% edges, and Taubin smoothing moves the vertex along the graph Laplacian.
% The positive lamda step removes local noise, and the negative mu step
% counteracts the shrinkage that ordinary Laplacian smoothing would cause.
%
% The smoothing is constrained by feat, ntype, and geom. Boundary passes
% smooth only boundary-surface vertices. Triple-line passes smooth only
% true non-manifold triple-line edges. Type 4 and type 14 quadruple junction
% vertices are used as geometric anchors and remain fixed during every pass.
%
% Usage:
%   Vout = graph_smooth_taubin(Vin,F,feat,ntype,geom,lamda,mu,num_iters);
%
% Inputs:
%   Vin       - [N x 3] input vertex coordinates.
%   F         - [M x 3] triangle connectivity using 1-based vertex indices.
%   feat      - Target feature to smooth: 'ext_triple', 'int_triple',
%               'bound', or another value for a general graph pass.
%   ntype     - [N x 1] vertex classification array. Types 3 and 13 are
%               triple-line vertices. Types 4 and 14 are quadruple junctions.
%   geom      - [N x 3] geometric per-coordinate smoothing weights. A row
%               with all zero values pins that vertex.
%   lamda     - Positive step size for the forward smoothing step.
%   mu        - Negative step size for the reverse anti-shrinkage step.
%   num_iters - Number of Taubin smoothing iterations.
%
% Outputs:
%   Vout      - [N x 3] smoothed vertex coordinates.
%
% Notes:
%   In triple-line passes, quadruple junctions are included in the feature
%   mask only so that real triple-line edges touching them can be detected.
%   They are then removed from the active set and are never moved.
%
% Change log:
%   2026/06 - Improve function graph_smooth_taubin.
%     Problem with the pristine function graph_smooth_taubin:
%       The original triple-line graph accepted any edge whose endpoints
%       had triple or quadruple labels. This included ordinary two-face
%       surface edges near quadruple points, allowed type 4/type 14
%       movement, and could over-pull nearby triple-line vertices.
%     
%     Corrections in the revised function graph_smooth_taubin:
%       Triple-line passes now keep only non-manifold edges shared by
%       more than two faces, use uniform connection weights, and pin
%       type 4/type 14 vertices during all smoothing passes.
%
%   2026/02 - Implement the surface smoothing method of XtalMesh in Matlab.
%             XtalMesh used Laplacian smoothing, but I prefer Taubin 
%             smoothing (without volume shrinkage). 
%
% Copyright (C) 2019-2026 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
%
% Project website: https://github.com/mjx888/im2mesh
%

    % Count mesh vertices once so later masks all use the same length.
    N = size(Vin, 1);

    % Feature passes need a stricter edge graph than regular surface passes.
    is_feature_pass = false;

    % ---------------------------------------------------------------------
    % Select vertices that may move during this smoothing pass.
    switch feat
        case 'ext_triple'
            is_feature_pass = true;

            % Exterior triple lines can terminate at exterior quadruple
            % junctions. Include type 14 for edge detection, but only allow
            % type 13 vertices to move.
            feature_mask = ismember(ntype, [13, 14]);
            active_mask = (ntype == 13);
        case 'int_triple'
            is_feature_pass = true;

            % Interior triple lines follow the same rule: type 4 vertices
            % help identify the line graph, but only type 3 vertices move.
            feature_mask = ismember(ntype, [3, 4]);
            active_mask = (ntype == 3);
        case 'bound'

            % Boundary smoothing is restricted to ordinary boundary-surface
            % vertices. Triple and quadruple features are handled separately.
            active_mask = ismember(ntype, [2, 12]);
        otherwise

            % A non-special pass may start from all vertices, then the global
            % anchor rules below remove fixed or geometrically blocked ones.
            active_mask = true(N, 1);
    end

    % Quadruple junctions are global anchors and never move in any pass.
    % Vertices with zero geometric weights are also pinned by constraint.
    is_quadruple = ismember(ntype, [4, 14]);
    active_mask = active_mask & ~is_quadruple & any(geom > 0, 2);
    active_idx = find(active_mask);
    N_act = length(active_idx);

    % If no vertices are allowed to move, return the input mesh unchanged.
    if N_act == 0, Vout = Vin; return; end

    % ---------------------------------------------------------------------
    % Construct adjacency for either a feature line or the surface mesh.
    if is_feature_pass
        % Build all undirected triangle edges. Counting repeated edges tells
        % how many surface triangles are incident to each mesh edge.
        edges = sort([
            F(:, [1, 2])
            F(:, [2, 3])
            F(:, [3, 1])
            ], 2);
        [unique_edges, ~, edge_index] = unique(edges, 'rows');
        edge_face_count = accumarray(edge_index, 1);

        % A true non-manifold triple-line edge is shared by more than two
        % triangles and has feature-type endpoints. This rejects accidental
        % two-face edges between nearby feature-labeled vertices.
        keep_edge = (edge_face_count > 2) & ...
            feature_mask(unique_edges(:, 1)) & ...
            feature_mask(unique_edges(:, 2));
        feature_edges = unique_edges(keep_edge, :);

        % Without a valid feature-line graph, this pass has nothing to smooth.
        if isempty(feature_edges), Vout = Vin; return; end

        % Store the feature-line graph in both directions so each edge acts as
        % an undirected graph connection during smoothing.
        i = [feature_edges(:, 1); feature_edges(:, 2)];
        j = [feature_edges(:, 2); feature_edges(:, 1)];
    else
        % For surface smoothing, keep faces touching at least one active
        % vertex. Fixed vertices on those faces remain available as neighbors.
        face_has_active = active_mask(F(:, 1)) | ...
            active_mask(F(:, 2)) | active_mask(F(:, 3));
        F_sub = F(face_has_active, :);

        % No incident faces means no graph neighbors for this pass.
        if isempty(F_sub), Vout = Vin; return; end

        % Convert each triangle to its three undirected graph edges, written
        % explicitly in both directions.
        i = [F_sub(:, 1); F_sub(:, 2); F_sub(:, 3); ...
             F_sub(:, 2); F_sub(:, 3); F_sub(:, 1)];
        j = [F_sub(:, 2); F_sub(:, 3); F_sub(:, 1); ...
             F_sub(:, 1); F_sub(:, 2); F_sub(:, 3)];
    end

    % Keep only graph rows for active vertices. Neighbor columns may still
    % refer to fixed vertices, because fixed vertices influence smoothing but
    % are not themselves updated.
    keep = active_mask(i);
    i_keep = i(keep);
    j_keep = j(keep);

    % Map global vertex ids to compact row ids for the active-vertex system.
    % Matrix columns remain in global coordinates to multiply the full V array.
    map_active = zeros(N, 1);
    map_active(active_idx) = 1:N_act;
    row_idx = map_active(i_keep);

    % Sparse adjacency from active vertices to all neighboring vertices.
    % spones removes duplicate triangle-edge entries and gives every graph
    % connection uniform weight, avoiding extra pull from quadruple vertices.
    A_act = sparse(row_idx, j_keep, 1, N_act, N);
    A_act = spones(A_act);

    % Degree and normalization for the unnormalized graph Laplacian:
    %     L(V) = degree .* V(active,:) - A * V
    % The 2*degree denominator preserves the legacy step scale and prevents
    % high-valence vertices from moving farther solely due to more neighbors.
    d_act = sum(A_act, 2);
    row_sum_abs_act = 2 * d_act;
    row_sum_abs_act(row_sum_abs_act == 0) = 1e-8;
    W_norm_act = 1 ./ row_sum_abs_act;

    % Combine the global Taubin parameters with per-coordinate geometric
    % weights. Zero geom entries suppress movement in their coordinate.
    pre_lamda_act = lamda * (geom(active_idx, :) .* W_norm_act);
    pre_mu_act    = mu    * (geom(active_idx, :) .* W_norm_act);

    % V contains all vertices because fixed neighbors still participate in
    % Laplacian evaluation. V_act caches only the rows that can be updated.
    V = Vin;
    V_act = V(active_idx, :);

    for it = 1:num_iters
        % Forward Laplacian step: move active vertices toward the local graph
        % average to remove jagged high-frequency surface noise.
        LV1_act = d_act .* V_act - A_act * V;
        V_half_act = V_act - pre_lamda_act .* LV1_act;
        V(active_idx, :) = V_half_act;

        % Reverse Taubin step: use the updated half-step mesh and a negative
        % mu value to compensate for shrinkage from the forward step.
        LV2_act = d_act .* V_half_act - A_act * V;
        V_new_act = V_half_act - pre_mu_act .* LV2_act;

        % Commit only active rows. Pinned vertices keep their input positions.
        V(active_idx, :) = V_new_act;
        V_act = V_new_act;
    end

    Vout = V;
end
