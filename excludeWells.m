% excludeWells.m - Exclude certain WellName.mat files from calculation

% Check if xlsx2mat.m has been run before
if (exist('expInfo.mat', 'file') ~= 2)
	if (exist('CLI', 'var') == 1)
		if (CLI == true(1))
			fprintf('expInfo.mat not found. Please run step 1 first.\n');
			return;
		end
	end
	fprintf('expInfo.mat not found. Please run xlsx2mat.m first.\n');
	return;
end

load('expInfo.mat');

useExcel = false;
while (1)
	useExcel = input('Use Excel file content? (Y/N): ', 's');
	% Detect if the length of the answer is 1, and if the answer is Y(y) or
	% N(n). If it is not, keep asking
	if (length(useExcel) == 1)
		if (useExcel(1) == 'y' || useExcel(1) == 'Y')
			useExcel = true(1);
			while (1)
				xlsxFilename = input('Input filename: ', 's');
				if (exist(xlsxFilename, 'file') ~= 2)
					xlsxFilename = sprintf('%s.xlsx', xlsxFilename);
					if (exist(xlsxFilename, 'file') == 2)
						break;
					end
					fprintf('File not found!\n');
				else
					break;
				end
			end
			break;
		elseif (useExcel(1) == 'n' || useExcel(1) == 'N')
			useExcel = false(1);
			break;
		end
	end
end

while (1)
	ifRmSysShift = input('Remove wells with significant systematic shift? (Y/N): ', 's');
	% Detect if the length of the answer is 1, and if the answer is Y(y) or
	% N(n). If it is not, keep asking
	if (length(ifRmSysShift) == 1)
		if (ifRmSysShift(1) == 'y' || ifRmSysShift(1) == 'Y')
			ifRmSysShift = true(1);
			while (1)
				nRmSysShift = input('How many frames of systematic shift is allowed: ');
				if (isnumeric(nRmSysShift))
					% If number of well is smaller than 1 or not integer, repeat question.
					if ((nRmSysShift >= 1) && (floor(nRmSysShift) == nRmSysShift))
						break;
					end
				end
			end
			if (QCdone)
				while (1)
					ifQC = input('Do you want to use QC data? (Y/N): ', 's');
					% Detect if the length of the answer is 1, and if the answer is Y(y) or
					% N(n). If it is not, keep asking
					if (length(ifQC) == 1)
						if (ifQC(1) == 'y' || ifQC(1) == 'Y')
							ifQC = true(1);
							break;
						elseif (ifQC(1) == 'n' || ifQC(1) == 'N')
							ifQC = false(1);
							break;
						end
					end
				end
			else
				ifQC = false(1);
			end
			break;
		elseif (ifRmSysShift(1) == 'n' || ifRmSysShift(1) == 'N')
			ifRmSysShift = false(1);
			break;
		end
	end
end

while (1)
	ifRmEmpty = input('Remove empty wells? (Y/N): ', 's');
	% Detect if the length of the answer is 1, and if the answer is Y(y) or
	% N(n). If it is not, keep asking
	if (length(ifRmEmpty) == 1)
		if (ifRmEmpty(1) == 'y' || ifRmEmpty(1) == 'Y')
			ifRmEmpty = true(1);
			break;
		elseif (ifRmEmpty(1) == 'n' || ifRmEmpty(1) == 'N')
			ifRmEmpty = false(1);
			break;
		end
	end
end

rmv = false(1, nWellProc);
if (useExcel)
	fprintf('Processing...');
	[~, ~, exclData] = xlsread(xlsxFilename, 'A1:X16');
	exclData(cellfun(@(x) all(isnan(x)), exclData)) = cellstr('');
	for i = 1:nWellProc
		row = char(rowList(i)) - 'A' + 1;
		col = str2double(char(colList(i)));
		if ((size(exclData, 1) - row >= 0) && (size(exclData, 2) - col >= 0))
			excl = char(exclData(row, col));
			if (~isempty(excl))
				rmv(i) = true;
			end
		end
	end
else
	fprintf('Please enter well label to be removed. One well per line.\n');
	fprintf('Format: 3-letter well label. E.g. A01, C12, F18\n');
	fprintf('Input ''end'' to end input.\n')
	while (1)
		wellName = input('', 's');
		if (length(wellName) == 3)
			if (strcmp(wellName, 'end'))
				break;
			end
			rowName = wellName(1);
			colName = wellName(2:3);
			temp = strcmp(rowList, rowName) && strcmp(colList, colName);
			rmv(temp) = true;
		end
	end
end

% Remove wells with significant systematic shift
for i = 1:nWellProc
	row = char(rowList(i)) - 'A' + 1;
	col = str2double(char(colList(i)));
	fileName = sprintf('Well%s%s.mat', char(rowList(i)), char(colList(i)));
	load(fileName);
	if (cellNum == 0)
		if (ifRmEmpty)
			rmv(i) = true;
		end
		continue;
	end
	if (ifRmSysShift)
		if (ifQC)
			dIntTraj = intTrajQC - wshift('2d', intTrajQC, [0 -1]);
			for k = 2:trjDuration
				sigmaDInt(k) = std(dIntTraj(:, k));
				muDInt(k) = mean(dIntTraj(:, k));
			end
			if (sum(sigmaDInt./abs(muDInt)<=1) > nRmSysShift)
				rmv(i) = true;
			end
		else
			dIntTraj = intTraj - wshift('2d', intTraj, [0 -1]);
			for k = 2:trjDuration
				sigmaDInt(k) = std(dIntTraj(:, k));
				muDInt(k) = mean(dIntTraj(:, k));
			end
			if (sum(sigmaDInt./abs(muDInt)<=1) > nRmSysShift)
				rmv(i) = true;
			end
		end
	end
end

rowList(rmv) = [];
colList(rmv) = [];
nWellProc = nWellProc - sum(rmv);

if (exist('tagList', 'var') == 1)
	tagList(rmv) = [];
	nTag = 1;
	uniqTag = tagList(1);
	for i = 2:nWellProc
		ifNewTag = 1;
		for j = 1:nTag
			if (strcmp(tagList(i), uniqTag(j)))
				ifNewTag = 0;
				break;
			end
		end
		if (ifNewTag)
			nTag = nTag + 1;
			uniqTag(nTag) = tagList(i);
		end
	end

	tagCount = zeros(1, nTag);
	for i = 1:nTag
		tagCount(i) = sum(strcmp(tagList, uniqTag(i)));
	end

	save('expInfo.mat', 'tagList', 'uniqTag', 'nTag', 'tagCount', 'rowList',...
		'colList', 'nWellProc', '-append');
else
	save('expInfo.mat', 'rowList', 'colList', 'nWellProc', '-append');
end

if (ifRmSysShift)
	save('expInfo.mat', 'ifRmSysShift', 'nRmSysShift', '-append');
end

if (useExcel)
	fprintf('Done!\n');
end