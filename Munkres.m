
%[tracks, newTrack_num, tracksRemoved_num] = hungarianAlgorithm(measurements, existingTracks)

%TODO:
%RESOLVE RECTANGULAR PROBLEMS, M MATRIX IS RECTANGULAR ()THERE ARE MORE MEASUREMENTS THAN TRACKS OR VICEVERSA


%hungarian algorithm
%rows-> measurements
%columns-> existing tracks
%not two measurements can be assigned to the same track

%measurements.range
%measurements.azimuth
%measurements.elevation
%measurements.velocity

%Initialize cost matrix
%n =size(existingTracks)
%M = zeros(n)
n=3;
%populate cost matrix to unit testing
%M=[40 60 15; 25 30 45; 55 30 25];
M = [15 15 0; 0 0 10; 5 5 0]; 

%step 1: row reduction: find minimum of each row
M_rowReduced = M;
for i = 1:n  %go through the rows
    minValueRow = min(M(i,:)); %minimum of that row
    for j = 1:n %go though the columns
        M_rowReduced(i,j) = M_rowReduced(i,j) - minValueRow ; 
    end
end

%step 2: column reduction: find minimum of each column
M_reduced = zeros(n);
for ii = 1:n  %go through the columns
    minValueColumn = min(M_rowReduced(:, ii)); %minimum of that column
    for jj = 1:n %go though the rows
        M_reduced(jj,ii) = M_rowReduced(jj,ii) - minValueColumn; 
    end
end

%% step3: count lines to cover all zeros
zeros_location = (M_reduced==0);
zeros_elements = find(zeros_location);
zeros_num = size(zeros_elements, 1);
zeros_matrix = zeros(2, zeros_num); % fila 1: fila del cero, fila 2: columna del cero, hay tantas columnas como ceros haya en la matriz
for k =1:zeros_num
    row_index = mod(zeros_elements(k),n); %it should be row index
    if (row_index == 0)
        row_index = n;
    end
    multiple = floor((zeros_elements(k)-1)/3); %it has to be 0 based, that´s why 1 is substracted
    column_index = multiple + 1; %matlab is not 0 based (1 based) 
    
        zeros_matrix(2,k) = column_index; %columna
        zeros_matrix(1,k) = row_index; %fila      
end

M_reduced_copy = M_reduced;
row_pp = 1;
rowsRemoved_num = 0;
for (pp= 1:n)
    [mode_value, o] = mode(M_reduced_copy(row_pp, :));  %o: occurrence
    if (mode_value == 0 && o > 1)
       M_reduced_copy(row_pp, :) = []; % remove row pp
       rowsRemoved_num = rowsRemoved_num + 1;
       rowsRemoved_record(rowsRemoved_num) = pp; %original row number removed
    else
        row_pp = row_pp + 1;
    end
    if (row_pp + rowsRemoved_num >= n+1)
        break
    end
end

column_pk = 1;
columnsRemoved_num = 0;
for (pk= 1:n)
    [mode_value, o] = mode(M_reduced_copy(:, column_pk));  %o: occurrence
    
    if ( (mode_value == 0 && o > 1) || (mode_value == 0 && o >= 1 && size(M_reduced_copy, 1) == 1) )
       M_reduced_copy(:, column_pk) = []; % remove row pp
       columnsRemoved_num = columnsRemoved_num + 1;
       columnsRemoved_record(columnsRemoved_num) = pk; %original column number removed
    else
        column_pk = column_pk + 1;
    end
    if (column_pk + columnsRemoved_num >= n+1)
        break
    end
end


%algorithm check -> M_reduced_copy shall not have any zero element
if (isempty(M_reduced_copy(M_reduced_copy == 0)) == 1);
    %everything is OK
else %exception
    disp("the lines counter to remove all zeros has failed")
end

lines_needed = rowsRemoved_num + columnsRemoved_num; % remove this

%% step 4 shift zeros

if (lines_needed == n)
    solutions_num = 1; %unique solution
    % do nothing
elseif (lines_needed < n)
    solutions_num = 2; % two or more 
    %shift zeros
    %to do
    %find minimum value uncovered by the remove rows and columns
    min_value = min(M_reduced_copy(:));
    %intersections_num = size(rowsRemoved_record, 2) * size(columnsRemoved_record, 2);
        for (ijk = 1: size(rowsRemoved_record, 2))
            for (kji = 1: size(columnsRemoved_record, 2))
                 M_reduced(rowsRemoved_record(ijk), columnsRemoved_record(kji)) = M_reduced(rowsRemoved_record(ijk), columnsRemoved_record(kji)) + 2*min_value; %twice the minimum  value is added because it will be substracted once later
            end
        end
        
    M_reduced((M_reduced ~= 0)) = M_reduced((M_reduced ~= 0)) - min_value;
    
    zeros_matrix = determineZeroRowColumn(M_reduced);
end

%step 5: lines to cover zeros = n
%choose zeros ensuring that each row and column only contains one chosen zero

rows = zeros_matrix(1,:);
columns = zeros_matrix(2,:);


%select zeros per row
rows_dummy = zeros_matrix(1,:);
f = 2;
while f > 1
    [currentMode, f] = mode(rows_dummy);
    if (f > 1)
        rows_dummy = rows_dummy(find(rows_dummy ~= currentMode));
    end
end
rowsWOMode = rows_dummy;
selectedZeros = zeros(2,n); %2-> row and column, n-> number of zeros selected
for kk=1:size(rowsWOMode,2)
 index = find(zeros_matrix(1,:)== rowsWOMode(kk));
 selectedZeros(1, kk) = zeros_matrix(1, index);
 selectedZeros(2, kk) = zeros_matrix(2, index);
end

%select zeros per column
columns_dummy = zeros_matrix(2,:);
f = 2;
while f > 1
    [currentMode_column, f] = mode(columns_dummy);
    if (f > 1)
        columns_dummy = columns_dummy(find(columns_dummy ~= currentMode_column));
    end
end
columnsWOMode = columns_dummy;
zerosDueToColumn = 0
for ij=1:size(columnsWOMode,2)
 index = find(zeros_matrix(2,:)== columnsWOMode(ij));
 element = [zeros_matrix(1, index); zeros_matrix(2, index)];
 differenceArray = element ~= selectedZeros;
 if (differenceArray == ones(2,n)) %check it is a genuine zero
    selectedZeros(1, kk + ij) = zeros_matrix(1, index);
    selectedZeros(2, kk + ij) = zeros_matrix(2, index);
    zerosDueToColumn = zerosDueToColumn + 1;
 end
end

zerosDueToRow = kk;
num_selectedZeros = zerosDueToRow + zerosDueToColumn,
if (zerosDueToRow + zerosDueToColumn == n)
    %do nothing, the n zeros has already been selected
else %we need to select another zero(s)
    full = [1:n]; 
    sample_rows = selectedZeros(1,:);
    idx = ismember(full, sample_rows);
    selectedZeros(1, num_selectedZeros+1:n) = full(~idx);
    
    %sample_columns = selectedZeros(2,:);
    %idx = ismember(full, sample_columns);
    %selectedZeros(2, n) = full(~idx);
    
    for ik = num_selectedZeros+1:n
        idx_ik = find(zeros_matrix(1,:) == selectedZeros(1, ik));
         %check zeros_matrix(2, idx_ik(1) is not already in selected zeros
         %if not choose zeros_matrix(2, idx_ik(2) etc 
        selectedZeros(2, ik) = zeros_matrix(2,idx_ik);
    end
    
end

%display, debug
M
M_reduced
selectedZeros