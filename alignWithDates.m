function alignedData = alignWithDates(dataTable, desiredDates)
    % Extract dates and data from the input table
    dataDates = dataTable.Var1; % Assume Var1 contains the dates
    dataValues = dataTable{:, 2:end}; % Extract the data columns

    % Initialize aligned data with NaN
    alignedData = nan(length(desiredDates), size(dataValues, 2));

    % Match dates in the input table with desiredDates
    [~, idxDesired, idxData] = intersect(desiredDates, dataDates, 'stable');

    % Fill in the aligned data for matched dates
    alignedData(idxDesired, :) = dataValues(idxData, :);
end