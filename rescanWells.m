% rescanWells.m - Rescan all .mat files

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

nWellProc = 0;
rowList = cell(1);
colList = cell(1);
for row = 65:80				% ASCII code for A ~ P (16 rows)
	for col = 1:24
		
		% Write filename
		if (col > 9)
			filename = sprintf('Well%c%d.mat', char(row), col);
		else
			filename = sprintf('Well%c0%d.mat', char(row), col);
		end
		
		% First see if the file exists. If not, directly go to next loop
		if (exist(filename, 'file') ~= 2)
			continue;
		end
		% Increment the number of wells processed.
		nWellProc = nWellProc + 1;
		
		% Store the name of the current well for future uses
		rowList(nWellProc) = cellstr(char(row));
		if (col > 9)
			colList(nWellProc) = cellstr(num2str(col));
		else
			colList(nWellProc) = cellstr(sprintf('0%d', col));
		end
	end
end

save('expInfo.mat', 'rowList', 'colList', 'nWellProc', '-append');

fprintf('Finished rescanning of %d .mat files!\n', nWellProc);
fprintf('Please re-run the tagging process to restore tags!\n');