function table = createLatexTable(data, rowNames, columnNames, caption, label, onlyData)

[rows, columns] = size(data);
if rows ~= length(rowNames)
    error('Rowlength and length of rowNames must match!');
end

if columns ~= length(columnNames)-1
    error('Number of columns must match number of columnNames');
end

prefix = '';
if ~onlyData
    prefix = ['\\begin{table}[ht]',...
                '\\centering'];


    str = '|';
    for i = 1: columns
       str = [str, 'c|']; 
    end

    prefix = [prefix, '\\begin{tabular}{', str, '}',...
              '\\hline'];


    tableColumnNames = '';
    for i = 1: columns
        columnName = columnNames(i);
        tableColumnNames = [tableColumnNames, strrep(columnName{1}),'\','\\'];
    end
end

tableData = '';
for row = 1:rows
    rowName = rowNames(row);
    tableData = [tableData, strrep(rowName{1}, '\', '\\')];
    for column = 1:columns
        value = data(row, column);
        if isa(value,'double')
            value = num2str(value);
        end
        tableData = [tableData, ' & ', value];
    end
    tableData = [tableData, '\\\\ \\hline\n'];
end
postfix = '';
if ~onlyData
    postfix = ['\\end{tabular}', '\\caption{', caption, '}', '\\label{tab:', label, '}', '\\end{table}'];
end

table = [prefix, tableData, postfix];
end