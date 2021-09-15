function [zeros_matrix] = determineZeroRowColumn(matrix)

n = size(matrix, 1);
zeros_location = (matrix==0);
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
