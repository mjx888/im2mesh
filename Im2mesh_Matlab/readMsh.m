function [nodes, elements] = readMshA(filename)
% readMshA: Reads an ASCII or binary Gmsh .msh file (Format 2.2).
% Extracts node coordinates and 4-node tetrahedral elements.
% Other element types, such as points, lines, and triangles, are ignored.
%
% Usage:
%   [nodes, elements] = readMshA('mesh.msh');
%
% Inputs:
%   filename - Path to an ASCII or binary Gmsh 2.2 .msh file
%
% Outputs:
%   nodes    - [N x 3] matrix containing x, y, z node coordinates
%   elements - [E x 4] matrix containing tetrahedron node IDs
%
% Notes:
%   1. The function automatically detects whether the .msh file is stored
%      in ASCII or binary format.
%   2. Only Gmsh file format version 2.2 is supported.
%   3. Binary files with little-endian or big-endian byte order are
%      supported.
%   4. The node IDs in elements are returned exactly as stored by Gmsh.
%      They are not renumbered to match the row numbers of nodes.
%
% Copyright (C) 2019-2026 by Jiexian Ma, mjx0799@gmail.com
% Distributed under the terms of the GNU General Public License (version 3)
%
% Project website: https://github.com/mjx888/im2mesh
%

    % ---------------------------------------------------------------------
    % Check the input filename.
    if ~(ischar(filename) || (isstring(filename) && isscalar(filename)))
        error('readMsh:InvalidFilename', ...
            'Filename must be a character vector or string scalar.');
    end

    % ---------------------------------------------------------------------
    % Open the beginning of the file in binary mode.
    % Binary mode can safely read the ASCII header of either file type.
    [fid, message] = fopen(filename, 'rb');
    if fid == -1
        error('readMsh:OpenFailed', 'Cannot open file: %s\n%s', ...
            char(filename), message);
    end
    % Close the file automatically if an error occurs.
    fileCleanup = onCleanup(@() fclose(fid));

    % ---------------------------------------------------------------------
    % Read the Gmsh mesh-format header.
    % Header format: version-number file-type data-size
    % file-type == 0: ASCII file
    % file-type == 1: binary file
    meshFormatSection = readLine(fid);
    formatValues = sscanf(readLine(fid), '%f');

    % Verify that the file is a supported Gmsh 2.2 mesh.
    if ~strcmp(meshFormatSection, '$MeshFormat')
        error('readMsh:InvalidFormat', ...
            'The file does not begin with a Gmsh $MeshFormat section.');
    end
    if numel(formatValues) ~= 3 || formatValues(1) ~= 2.2
        error('readMsh:UnsupportedVersion', ...
            'Only Gmsh file format version 2.2 is supported.');
    end

    fileType = formatValues(2);
    % Close the header stream before reopening the file in the proper mode.
    delete(fileCleanup);

    % ---------------------------------------------------------------------
    % Select the appropriate parser according to the file-type flag.
    switch fileType
        case 0
            [nodes, elements] = readAsciiFile(filename);
        case 1
            [nodes, elements] = readBinaryFile(filename);
        otherwise
            error('readMsh:InvalidFileType', ...
                'Invalid Gmsh file type value: %g.', fileType);
    end
end

function [nodes, elements] = readAsciiFile(filename)
% readAsciiFile: Read nodes and tetrahedral elements from an ASCII file.

    % ---------------------------------------------------------------------
    % Open the ASCII .msh file for reading.
    [fid, message] = fopen(filename, 'r');
    if fid == -1
        error('readMsh:OpenFailed', 'Cannot open file: %s\n%s', ...
            char(filename), message);
    end
    fileCleanup = onCleanup(@() fclose(fid));

    % Initialize outputs before searching for the Nodes and Elements
    % sections.
    nodes = [];
    elements = zeros(0, 4);

    % ---------------------------------------------------------------------
    % Read the file section by section.
    while ~feof(fid)
        section = readLine(fid);
        switch section
            case '$Nodes'
                % ASCII node format:
                % node-number x-coordinate y-coordinate z-coordinate
                numNodes = readPositiveIntegerLine(fid, 'node count');
                nodeData = fscanf(fid, '%f', [4, numNodes]).';
                if size(nodeData, 1) ~= numNodes
                    error('readMsh:UnexpectedEndOfFile', ...
                        'Unexpected end of file while reading ASCII nodes.');
                end
                % Discard node numbers and keep only x, y, z coordinates.
                nodes = nodeData(:, 2:4);

            case '$Elements'
                % Read all element records and retain type 4 elements.
                % Gmsh element type 4 is a 4-node tetrahedron.
                numElements = readPositiveIntegerLine(fid, 'element count');
                elements = readAsciiElements(fid, numElements);

            case '$ElementData'
                % Node and element connectivity have already been read.
                % Ignore optional result data that follows.
                break;
        end
    end

    % Check that both required mesh components were found.
    validateMesh(nodes, elements);
end

function elements = readAsciiElements(fid, numElements)
% readAsciiElements: Read ASCII element records and extract tetrahedra.

    % Preallocate for the maximum possible number of tetrahedra.
    elements = zeros(numElements, 4);
    tetCount = 0;

    % ---------------------------------------------------------------------
    % ASCII element format:
    % element-number element-type number-of-tags tags... node-numbers...
    for elementIndex = 1:numElements
        values = sscanf(readLine(fid), '%f');
        if numel(values) < 3
            error('readMsh:InvalidElement', ...
                'Invalid ASCII element record number %d.', elementIndex);
        end

        if values(2) == 4
            % Tags occur before the node-number list, so use the number of
            % tags to locate the four tetrahedron node IDs.
            numTags = values(3);
            nodesStart = 4 + numTags;
            if numel(values) < nodesStart + 3
                error('readMsh:InvalidElement', ...
                    'Invalid tetrahedral element record number %d.', ...
                    elementIndex);
            end
            tetCount = tetCount + 1;
            elements(tetCount, :) = values(nodesStart:nodesStart + 3).';
        end
    end

    % Remove unused rows created by preallocation.
    elements = elements(1:tetCount, :);
end

function [nodes, elements] = readBinaryFile(filename)
% readBinaryFile: Read nodes and tetrahedral elements from a binary file.

    % ---------------------------------------------------------------------
    % First open the file as little-endian. The endianness check integer in
    % the mesh-format section determines whether it must be reopened.
    [fid, message] = fopen(filename, 'rb', 'ieee-le');
    if fid == -1
        error('readMsh:OpenFailed', 'Cannot open file: %s\n%s', ...
            char(filename), message);
    end
    fileCleanup = onCleanup(@() fclose(fid));

    % Read the binary mesh header and its 4-byte endianness check integer.
    requireLine(fid, '$MeshFormat');
    readLine(fid);
    endianCheck = fread(fid, 1, 'int32=>int32');
    isBigEndian = false;

    % The check integer must equal 1. If its swapped value equals 1, the
    % file uses big-endian byte order and must be reopened accordingly.
    if endianCheck ~= 1
        if swapbytes(endianCheck) ~= 1
            error('readMsh:InvalidEndianCheck', ...
                'Invalid binary endianness check value.');
        end

        delete(fileCleanup);
        [fid, message] = fopen(filename, 'rb', 'ieee-be');
        if fid == -1
            error('readMsh:OpenFailed', 'Cannot reopen file: %s\n%s', ...
                char(filename), message);
        end
        fileCleanup = onCleanup(@() fclose(fid));
        isBigEndian = true;

        requireLine(fid, '$MeshFormat');
        readLine(fid);
        endianCheck = fread(fid, 1, 'int32=>int32');
        if endianCheck ~= 1
            error('readMsh:InvalidEndianCheck', ...
                'Invalid binary endianness check value.');
        end
    end

    % Finish reading the mesh-format section.
    requireNextNonemptyLine(fid, '$EndMeshFormat');

    % Initialize outputs before searching for binary mesh sections.
    nodes = [];
    elements = zeros(0, 4);

    % ---------------------------------------------------------------------
    % Read the file section by section. Text section markers surround the
    % binary node and element records in a Gmsh 2.2 file.
    while ~feof(fid)
        section = readLine(fid);
        if isempty(section)
            continue;
        end

        switch section
            case '$Nodes'
                % Read all binary node records in one operation.
                numNodes = readPositiveIntegerLine(fid, 'node count');
                nodes = readBinaryNodes(fid, numNodes, isBigEndian);
                requireNextNonemptyLine(fid, '$EndNodes');

            case '$Elements'
                % Binary elements are grouped into blocks by element type
                % and number of tags.
                numElements = readPositiveIntegerLine(fid, 'element count');
                elements = readBinaryElements(fid, numElements);
                requireNextNonemptyLine(fid, '$EndElements');

            case '$ElementData'
                % Ignore optional result data after mesh connectivity.
                break;

            otherwise
                % Skip optional text sections not needed by this function.
                if startsWith(section, '$') && ~startsWith(section, '$End')
                    skipSection(fid, section);
                end
        end
    end

    % Check that both required mesh components were found.
    validateMesh(nodes, elements);
end

function nodes = readBinaryNodes(fid, numNodes, isBigEndian)
% readBinaryNodes: Read the fixed-length binary node records.

    % Each node record contains:
    %   one 4-byte integer node ID + three 8-byte double coordinates.
    bytesPerNode = 4 + 3 * 8;
    raw = fread(fid, [bytesPerNode, numNodes], 'uint8=>uint8');
    if numel(raw) ~= bytesPerNode * numNodes
        error('readMsh:UnexpectedEndOfFile', ...
            'Unexpected end of file while reading binary nodes.');
    end

    % Discard the first four bytes of each record, which contain node IDs.
    % Convert the remaining bytes into x, y, z double values.
    coordinateBytes = reshape(raw(5:end, :), [], 1);
    coordinateValues = typecast(coordinateBytes, 'double');
    if isBigEndian
        % typecast uses the computer byte order, so swap big-endian data.
        coordinateValues = swapbytes(coordinateValues);
    end
    % Store one node coordinate triplet per matrix row.
    nodes = reshape(coordinateValues, 3, numNodes).';
end

function elements = readBinaryElements(fid, numElements)
% readBinaryElements: Read binary element blocks and extract tetrahedra.

    % Preallocate for the maximum possible number of tetrahedra.
    elements = zeros(numElements, 4);
    tetCount = 0;
    elementsRead = 0;

    % ---------------------------------------------------------------------
    % Each binary element block starts with three 4-byte integers:
    % element-type, number-of-elements-in-block, and number-of-tags.
    while elementsRead < numElements
        blockHeader = fread(fid, 3, 'int32=>double');
        if numel(blockHeader) ~= 3
            error('readMsh:UnexpectedEndOfFile', ...
                'Unexpected end of file while reading an element block.');
        end

        elementType = blockHeader(1);
        blockCount = blockHeader(2);
        numTags = blockHeader(3);
        numNodes = nodesPerElement(elementType);

        % Each element record contains one element ID, zero or more tags,
        % and a fixed number of node IDs determined by its element type.
        valuesPerElement = 1 + numTags + numNodes;
        block = fread(fid, [valuesPerElement, blockCount], 'int32=>double');

        if numel(block) ~= valuesPerElement * blockCount
            error('readMsh:UnexpectedEndOfFile', ...
                'Unexpected end of file while reading element type %d.', ...
                elementType);
        end

        if elementType == 4
            % Skip the element ID and tags, then copy the four node IDs.
            nodeRows = (2 + numTags):(5 + numTags);
            elements(tetCount + (1:blockCount), :) = block(nodeRows, :).';
            tetCount = tetCount + blockCount;
        end
        elementsRead = elementsRead + blockCount;
    end

    % Remove unused rows belonging to ignored element types.
    elements = elements(1:tetCount, :);
end

function numNodes = nodesPerElement(elementType)
% nodesPerElement: Return the node count for a Gmsh 2.2 element type.

    % Entry n stores the number of nodes used by Gmsh element type n.
    nodeCounts = [2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, ...
        1, 8, 20, 15, 13, 9, 10, 12, 15, 15, 21, 4, 5, 6, 20, 35, 56];
    if elementType < 1 || elementType > numel(nodeCounts)
        error('readMsh:UnsupportedElementType', ...
            'Unsupported Gmsh element type: %d.', elementType);
    end
    numNodes = nodeCounts(elementType);
end

function validateMesh(nodes, elements)
% validateMesh: Verify that nodes and tetrahedra were successfully read.

    if isempty(nodes)
        error('readMsh:MissingNodes', 'No nodes were found in the file.');
    end
    if isempty(elements)
        error('readMsh:MissingTetrahedra', ...
            'No 4-node tetrahedral elements were found in the file.');
    end
    fprintf('Loaded %d nodes and %d elements.\n', size(nodes, 1), ...
        size(elements, 1));
end

function count = readPositiveIntegerLine(fid, description)
% readPositiveIntegerLine: Read and validate a nonnegative integer line.

    count = sscanf(readLine(fid), '%d', 1);
    if isempty(count) || count < 0
        error('readMsh:InvalidCount', 'Invalid %s.', description);
    end
end

function requireLine(fid, expected)
% requireLine: Verify that the next line equals the expected text.

    actual = readLine(fid);
    if ~strcmp(actual, expected)
        error('readMsh:InvalidFormat', 'Expected "%s", found "%s".', ...
            expected, actual);
    end
end

function requireNextNonemptyLine(fid, expected)
% requireNextNonemptyLine: Verify the next nonempty section marker.
% Some Gmsh writers add a newline between binary data and the next marker,
% while others write the marker immediately after the binary data.

    actual = readLine(fid);
    while isempty(actual) && ~feof(fid)
        actual = readLine(fid);
    end
    if ~strcmp(actual, expected)
        error('readMsh:InvalidFormat', 'Expected "%s", found "%s".', ...
            expected, actual);
    end
end

function line = readLine(fid)
% readLine: Read one trimmed text line from the current file position.

    line = fgetl(fid);
    if ~ischar(line)
        line = '';
    else
        line = strtrim(line);
    end
end

function skipSection(fid, sectionStart)
% skipSection: Skip an optional text section not used by this function.

    % For example, sectionStart '$PhysicalNames' ends at
    % '$EndPhysicalNames'.
    sectionEnd = ['$End', sectionStart(2:end)];
    while ~feof(fid)
        if strcmp(readLine(fid), sectionEnd)
            return;
        end
    end
    error('readMsh:InvalidFormat', 'Missing section terminator "%s".', ...
        sectionEnd);
end
