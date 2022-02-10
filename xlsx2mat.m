% xlsx2mat.m - Reads Excel files from the macro output and outputs .mat files
% Assuming no gaps within a single trajectory
% Saves equivalent diameter, mean intensity, time point and actual time
% Calculates and saves number of cells, start and end point of each
% trajectories
% Dissects the raw data into trajectory matricies according to user
% requirement

% Number of wells. Input by user
while (1)
	nWell = input('Number of wells to process: ');
	
	% Detect if number of well is a number
	% If it is not, keep asking the number of wells
	if (isnumeric(nWell))
		% If number of well is smaller than 1 or not integer, repeat question.
		if (nWell < 1)
			fprintf('Number of wells has to be at least 1!\n');
		elseif (floor(nWell) ~= nWell)
			fprintf('Number of wells has to be an integer!\n');
		else
			break;
		end
	else
		fprintf('Number of wells has to be a number!\n');
	end
end

% At which frame should the dataset be starting at?
while(1)
	trjStart = input('Starting frame #: ');
	
	% Detect if starting frame # is a number
	% If it is not, keep asking the number of wells
	if (isnumeric(trjStart))
		% If starting frame # is smaller than 1 or not integer, repeat question.
		if (trjStart < 1)
			fprintf('Starting frame # has to be at least 1!\n');
		elseif (floor(trjStart) ~= trjStart)
			fprintf('Starting frame # has to be an integer!\n');
		else
			break;
		end
	else
		fprintf('Starting frame # has to be a number!\n');
	end
end

% How long of a span in each trajectory should be included in the dataset?
while(1)
	trjDuration = input('Number of frames to process: ');
	
	% Detect if number of frames is a number
	% If it is not, keep asking the number of wells
	if (isnumeric(trjDuration))
		% If number of frames is smaller than 1 or not integer, repeat question.
		if (trjDuration < 1)
			fprintf('Number of frames has to be at least 1!\n');
		elseif (floor(trjDuration) ~= trjDuration)
			fprintf('Number of frames has to be an integer!\n');
		else
			break;
		end
	else
		fprintf('Number of frames has to be a number!\n');
	end
end

% Error checking variable
err = false(1);
% Number of wells that has been processed
nWellProc = 0;

% Stores the columns and rows of all wells
rowList = cell(1, nWell);
colList = cell(1, nWell);

% Step through all wells of a 384-well plate
for row = 65:80				% ASCII code for A ~ P (16 rows)
	% Check if the given number of wells has been processed
	if (nWellProc == nWell)
		break;
	end
	for col = 1:24
		% Check if the given number of wells has been processed
		if (nWellProc == nWell)
			break;
		end
		
		% Write filename
		if (col > 9)
			filename = sprintf('Well%c%d_.xlsx', char(row), col);
		else
			filename = sprintf('Well%c0%d_.xlsx', char(row), col);
		end
		
		% First see if the file exists. If not, directly go to next loop
		if (exist(filename, 'file') ~= 2)
			continue;
		end
		fprintf('Reached well #%d/%d! Processing... ', nWellProc + 1, nWell);
		
		% Try to open the current well Excel sheet
		try
			[~, ~, raw] = xlsread(filename, 'Data');
		catch e
			if strcmp(e.identifier, 'MATLAB:xlsread:WorksheetNotFound')
				% If the worksheet 'Data' does not exist, there may be
				% potential problems with the Excel file. Abort the program
				fprintf('Cannot find sheet ''Data'' in file %s!\n', filename);
				err = true(1);
				break;
			else
				% If other errors occurred, abort the program as well
				fprintf('Unknown error occurred reading file %s!\n', filename);
				err = true(1);
				break;
			end
		end
		
		% If no error occurred during xlsread(), increment the number of
		% wells processed.
		nWellProc = nWellProc + 1;
		
		% Store the name of the current well for future uses
		rowList(nWellProc) = cellstr(char(row));
		if (col > 9)
			colList(nWellProc) = cellstr(num2str(col));
		else
			colList(nWellProc) = cellstr(sprintf('0%d', col));
		end

		% Determine columns with data of interest
		titleRow = raw(1, :);
		
		% Compare the clumn titles to find time, equivalent diameter, mean
		% intensity, and time point columns
		timeCol = find(strcmp(titleRow, 'Time [s]'));
		eqDiamCol = find(strcmp(titleRow, sprintf('EqDiameter [%sm]', 181))); % char(181) = [greek letter mu]
		intensityCol = find(strcmp(titleRow, 'MeanIntensity'));
		sumIntCol = find(strcmp(titleRow, 'SumIntensity'));
		timePtCol = find(strcmp(titleRow, 'ND.T'));
		multiCol = find(strcmp(titleRow, 'ND.M'));
		IDcol = find(strcmp(titleRow, 'ID'));
		
		% Extract numbers from the columns
		rawTime = raw(2:end, timeCol);
		eqDiam = cell2mat(raw(2:end, eqDiamCol))';				% in µm
		intensity = cell2mat(raw(2:end, intensityCol))';
		sumInt = cell2mat(raw(2:end, sumIntCol))';
		timePt = cell2mat(raw(2:end, timePtCol))';
		multiPoint = cell2mat(raw(2:end, multiCol))';
		ID = cell2mat(raw(2:end, IDcol))';
		
		% Detect time data type
		time = zeros(1, length(rawTime));
		for i = 1:length(rawTime)
			isString = false;
			tempNum = cell2mat(rawTime(i));
			if (isnumeric(tempNum) == 0)
				tempString = char(rawTime(i));
				ds = str2double(tempString(1));
				hrs = str2double(tempString(4:5));
				mins = str2double(tempString(8:9));
				time(i) = hrs * 3600 + mins * 60 + ds * 24 * 3600;
			else
				time(i) = tempNum * 24 * 60 * 60;
			end
		end

		% Separating single cell trajectories
		dMulti = wshift('1D', multiPoint, 1) - multiPoint;
		dID = wshift('1D', ID, 1) - ID;
		endOfTrajBool = 1 - ((dID == 0) & (dMulti == 0));
		endOfTrajIndex = find(endOfTrajBool);
		startOfTrajIndex = [1 endOfTrajIndex(1:end - 1) + 1];
		totalCellNum = length(endOfTrajIndex);
		
		% Calculating mean time interval between frames
		dt = wshift('1D', time, 1) - time;
		dt(endOfTrajBool == 1) = [];
		meanFrameDuration = mean(dt);
		
		% Process the trajectories to only include those that have the 
		% required length
		cellNum = 0;
		intTraj = zeros(totalCellNum, trjDuration);
		sumIntTraj = intTraj;
		eqDiamTraj = intTraj;
		timeTraj = intTraj;
		timePtTraj = intTraj;
		for i = 1:totalCellNum
			if (timePt(startOfTrajIndex(i)) > trjStart)
				continue;
			end
			if (timePt(endOfTrajIndex(i)) < (trjStart + trjDuration - 1))
				continue;
			end
			cellNum = cellNum + 1;
			tempTimePt = timePt(startOfTrajIndex(i):endOfTrajIndex(i));
			tempIntensity = intensity(startOfTrajIndex(i):endOfTrajIndex(i));
			tempSumInt = sumInt(startOfTrajIndex(i):endOfTrajIndex(i));
			tempEqDiam = eqDiam(startOfTrajIndex(i):endOfTrajIndex(i));
			tempTime = time(startOfTrajIndex(i):endOfTrajIndex(i));
			intTraj(cellNum, :) = tempIntensity((tempTimePt >= trjStart) ...
				& (tempTimePt < trjStart + trjDuration));
			sumIntTraj(cellNum, :) = tempSumInt((tempTimePt >= trjStart) ...
				& (tempTimePt < trjStart + trjDuration));
			eqDiamTraj(cellNum, :) = tempEqDiam((tempTimePt >= trjStart) ...
				& (tempTimePt < trjStart + trjDuration));
			timeTraj(cellNum, :) = tempTime((tempTimePt >= trjStart) ...
				& (tempTimePt < trjStart + trjDuration));
			timePtTraj(cellNum, :) = tempTimePt((tempTimePt >= trjStart) ...
				& (tempTimePt < trjStart + trjDuration));
		end
		intTraj(cellNum + 1:end, :) = [];
		sumIntTraj(cellNum + 1:end, :) = [];
		eqDiamTraj(cellNum + 1:end, :) = [];
		timeTraj(cellNum + 1:end, :) = [];
		timePtTraj(cellNum + 1:end, :) = [];
		
		% Count cell number at each time frame
		for i = max(timePt):-1:min(timePt)
			nCell(i) = sum(timePt == i);
		end
		
		% Save to a .mat file
		if (col > 9)
			filename = sprintf('Well%c%d.mat', char(row), col);
		else
			filename = sprintf('Well%c0%d.mat', char(row), col);
		end
		save(filename, 'time', 'eqDiam', 'intensity', 'timePt', ...
			'endOfTrajIndex', 'startOfTrajIndex', 'totalCellNum', ...
			'trjStart', 'trjDuration', 'intTraj', 'eqDiamTraj', 'timeTraj', ...
			'timePtTraj', 'cellNum', 'meanFrameDuration', 'nCell', 'sumIntTraj');
		fprintf('Done!\n');
	end
	
	% If encountered error, abort the program
	if (err)
		break;
	end
end

if (~err)
	if (nWellProc == nWell)
		fprintf('Completed output of %d .mat file!\n', nWellProc);
	elseif (nWellProc < nWell)
		fprintf('Only %d .mat file was created before reaching the end of the plate!\n', ...
			nWellProc);
	else
		fprintf('%d (> %d) .mat file was created. This should not happen!\n', ...
			nWellProc, nWell);
	end
	
	% Flag the data for not yet passed QC
	QCdone = false(1);
	
	% Flag the data for not being merged according to tags
	merged = false;
	save('expInfo.mat', 'trjStart', 'trjDuration', 'nWell', 'nWellProc', ...
		'rowList', 'colList', 'QCdone', 'merged');
else
	fprintf('Aborting...\n');
end
