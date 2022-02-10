% tagWells.m - Add tags to WellName.mat files

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
if (exist('tag.xlsx', 'file') == 2)
	while (1)
		useExcel = input('Use tag.xlsx content for tags? (Y/N): ', 's');
		% Detect if the length of the answer is 1, and if the answer is Y(y) or
		% N(n). If it is not, keep asking
		if (length(useExcel) == 1)
			if (useExcel(1) == 'y' || useExcel(1) == 'Y')
				useExcel = true(1);
				fprintf('Processing...');
				[~, ~, tagData] = xlsread('tag.xlsx');
				break;
			elseif (useExcel(1) == 'n' || useExcel(1) == 'N')
				useExcel = false(1);
				break;
			end
		end
	end
end

tagList = cell(1, nWellProc);
for i = 1:nWellProc
	% Write filename
	filename = sprintf('Well%s%s.mat', char(rowList(i)), char(colList(i)));
	
	% Write well name
	wellName = sprintf('Well%s%s', char(rowList(i)), char(colList(i)));
	
	% Load data
	load(filename);
	
	if (useExcel)
		% Read tag and store
		row = char(rowList(i)) - 'A' + 1;
		col = str2double(char(colList(i)));
		try
			tag = char(tagData(row, col));
		catch
			tag = cell2mat(tagData(row, col));
			if isnan(tag)
				tag = '';
			else
				tag = num2str(tag);
			end
		end
	else
		% Ask for tag and store
		tag = input(sprintf('Tag for well %s: ', wellName), 's');
	end
	save(filename, 'tag', '-append');
	tagList(i) = cellstr(tag);
end

% Detect how many unique tags exist
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

if (useExcel)
	fprintf(' Done!\n');
end

% Save tag list into expInfo.mat
save('expInfo.mat', 'tagList', 'uniqTag', 'nTag', 'tagCount', '-append');