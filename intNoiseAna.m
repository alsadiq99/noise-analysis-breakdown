% intNoiseAna.m - Noise analysis for intensity trajectory data
% Currently only works for isoclonal data

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

% See if noise has been processed. If not, create an empty structure. If
% so, load old file and append to it.
if (exist('intNoiseInfo', 'var') ~= 1)
	intNoiseInfo = struct();
	intNoiseInfo.existQC = false(1);
	intNoiseInfo.existNonQC = false(1);
elseif (isstruct(intNoiseInfo) == 0)
	intNoiseInfo = struct();
	intNoiseInfo.existQC = false(1);
	intNoiseInfo.existNonQC = false(1);
end

% Check if merged tagged data exist
if (merged)
	while (1)
		useMerged = input('Use merged dataset according to tags? (Y/N): ', 's');
		% Detect if the length of the answer is 1, and if the answer is Y(y) or
		% N(n). If it is not, keep asking
		if (length(useMerged) == 1)
			if (useMerged(1) == 'y' || useMerged(1) == 'Y')
				useMerged = true(1);
				break;
			elseif (useMerged(1) == 'n' || useMerged(1) == 'N')
				useMerged = false(1);
				break;
			end
		end
	end
else
	useMerged = false;
end

% Ask the user if the data is polyclonal
while (1)
	isPoly = input('Is the data from polyclone cells? (Y/N): ', 's');
	% Detect if the length of the answer is 1, and if the answer is Y(y) or
	% N(n). If it is not, keep asking
	if (length(isPoly) == 1)
		if (isPoly(1) == 'y' || isPoly(1) == 'Y')
			isPoly = true(1);
			% Ask for number of clusters
			while (1)
				nCluster = input('Number of clusters: ');
				% Detect if answer is a number, is at least one, and is an integer
				% If it is not, keep asking
				if (isnumeric(nCluster) && nCluster >= 1 && floor(nCluster) == nCluster)
					break;
				end
			end
			break;
		elseif (isPoly(1) == 'n' || isPoly(1) == 'N')
			isPoly = false(1);
			break;
		end
	end
end

% Detects if QC was done, and asks user whether to use QC data or not
if (QCdone)
	while (1)
		useQC = input('Use QC data? (Y/N): ', 's');
		% Detect if the length of the answer is 1, and if the answer is Y(y) or
		% N(n). If it is not, keep asking
		if (length(useQC) == 1)
			if (useQC(1) == 'y' || useQC(1) == 'Y')
				useQC = true(1);
				break;
			elseif (useQC(1) == 'n' || useQC(1) == 'N')
				useQC = false(1);
				break;
			end
		end
	end
else
	useQC = false(1);
end

if (useQC)
	intNoiseInfo.existQC = true;
else
	intNoiseInfo.existNonQC = true;
end
if (isPoly)
	intNoiseInfo.isPoly = true;
	intNoiseInfo.nCluster = nCluster;
else
	intNoiseInfo.isPoly = false;
end

j = 1;
while (1)
	% Load data
	if (useMerged)
		% Check if all tags are finished
		if (j > nTag)
			break;
		end
		
		% Write filename
		filename = sprintf('%s.mat', char(uniqTag(j)));
		
		% Write well name
		wellName = char(uniqTag(j));
		
		% Load data
		load(filename);
		fprintf('Reached tag #%d/%d. Processing... ', j, nTag);
	else
		% Check if all wells are finished
		if (j > nWellProc)
			break;
		end
		
		% Write filename
		filename = sprintf('Well%s%s.mat', char(rowList(j)), char(colList(j)));
		
		% Write well name
		wellName = sprintf('Well%s%s', char(rowList(j)), char(colList(j)));
		
		% Load data
		load(filename);
		fprintf('Reached well #%d/%d. Processing... ', j, nWellProc);
	end
	j = j + 1;
	
	if (cellNum == 0)
		fprintf('No cells detected!\n');
		continue;
	end
	
	if (isPoly)
		% Divide up the clusters using end point clustering
		if (useQC)
			intTrajQCEnd = intTrajQC(:, end);
		else
			intTrajQCEnd = intTraj(:, end);
		end
		[~, tempInd] = sort(intTrajQCEnd);
		if (useQC)
			intClustIndQC = zeros(1, cellNumQC);
			for i = 1:cellNumQC
				intClustIndQC(tempInd(i)) = ceil(i / cellNumQC * nCluster);
			end
		else
			intClustIndQC = zeros(1, cellNum);
			for i = 1:cellNum
				intClustIndQC(tempInd(i)) = ceil(i / cellNum * nCluster);
			end
		end
		
		% Noise processing
		for i = 1:nCluster
			if (useQC)
				tempTraj = intTrajQC(intClustIndQC == i, :);
			else
				tempTraj = intTraj(intClustIndQC == i, :);
			end
			nTraj = sum(intClustIndQC == i);
			
			% 
			meanIntTrend = mean(tempTraj, 1);
			intTrajQCdet = tempTraj - repmat(meanIntTrend, [nTraj 1]);
			meanInt = mean(intTrajQCdet, 2);
			intNoiseQC = intTrajQCdet - repmat(meanInt, [1 trjDuration]);
			
			% Calculate mean population intensity
			intMeanQC = mean(mean(tempTraj));
			
			% Calculate zero crossing
			temp = intNoiseQC(1:end-1) .* wshift('1D', intNoiseQC(1:end-1), 1);
			intZCrossQC = sum(temp <= 0) / trjDuration / meanFrameDuration * 3600;
			
			% Calculate the ACF for intensity noise with QC
			totalAcf = zeros(1, trjDuration);
			intAcfTrajQC = zeros(nTraj, trjDuration);
			for k = 1:nTraj
				[acf, lag] = xcorr(intNoiseQC(k, :), 'biased');
				acf(lag < 0) = [];
				totalAcf = totalAcf + acf;
				intAcfTrajQC(k, :) = acf;
			end
			totalAcf = totalAcf / nTraj;
			
			% Variance (sigma2) is the zero-lag value of ACF
			intSigma2QC = totalAcf(1);
			
			% Cross correlation distance. Use squareform() to convert into
			% square matrix
			intNoiseXCorrQC = pdist(intNoiseQC, @xcorrDist);
			
			% Interpolate ACF linearly and find half maximum value
			try
				intTauHalfQC = interp1(totalAcf, 1:trjDuration, 0.5 * max(totalAcf));
				intTauHalfQC = intTauHalfQC * meanFrameDuration / 3600;
				intCV2QC = intSigma2QC / intMeanQC ^ 2;
			catch
				intTauHalfQC = -1;
				intCV2QC = -1;
			end
			
			% Create structure to hold information for current well
			% Check if this well has been processed for noise. If not,
			% create a new field for the well.
			if (~isfield(intNoiseInfo, wellName))
				intNoiseInfo.(wellName) = struct();
			end
			
			% Create struct array for noise data
			if (isfield(intNoiseInfo.(wellName), 'polyClust'))
				if (isstruct(intNoiseInfo.(wellName).polyClust))
					if (length(intNoiseInfo.(wellName).polyClust) ~= nCluster)
						intNoiseInfo.(wellName).polyClust = repmat(struct(), [1, nCluster]);
					end
				else
					intNoiseInfo.(wellName).polyClust = repmat(struct(), [1, nCluster]);
				end
			else
				intNoiseInfo.(wellName).polyClust = repmat(struct(), [1, nCluster]);
			end
			
			if (useQC)
				intNoiseInfo.(wellName).polyClust(i).intMeanQC = intMeanQC;
				intNoiseInfo.(wellName).polyClust(i).intSigma2QC = intSigma2QC;
				intNoiseInfo.(wellName).polyClust(i).intTauHalfQC = intTauHalfQC;
				intNoiseInfo.(wellName).polyClust(i).intCV2QC = intCV2QC;
				intNoiseInfo.(wellName).polyClust(i).intAcfTrajQC = intAcfTrajQC;
				intNoiseInfo.(wellName).polyClust(i).intNoiseTrajQC = intNoiseQC;
				intNoiseInfo.(wellName).polyClust(i).intZCrossQC = intZCrossQC;
				intNoiseInfo.(wellName).polyClust(i).intNoiseXCorrQC = intNoiseXCorrQC;
			else
				intNoiseInfo.(wellName).polyClust(i).intMean = intMeanQC;
				intNoiseInfo.(wellName).polyClust(i).intSigma2 = intSigma2QC;
				intNoiseInfo.(wellName).polyClust(i).intTauHalf = intTauHalfQC;
				intNoiseInfo.(wellName).polyClust(i).intCV2 = intCV2QC;
				intNoiseInfo.(wellName).polyClust(i).intAcfTraj = intAcfTrajQC;
				intNoiseInfo.(wellName).polyClust(i).intNoiseTraj = intNoiseQC;
				intNoiseInfo.(wellName).polyClust(i).intZCross = intZCrossQC;
				intNoiseInfo.(wellName).polyClust(i).intNoiseXCorr = intNoiseXCorrQC;
			end
		end
		% Save noise informatino to .mat file
		polyClust = intNoiseInfo.(wellName).polyClust;
		if (useQC)
			save(filename, 'polyClust', 'intClustIndQC', '-append');
		else
			intClustInd = intClustIndQC;
			save(filename, 'polyClust', 'intClustInd', '-append');
		end
		fprintf('Done!\n');
	else
		% Calculate noise using QC data
		if (useQC)
			meanIntTrend = mean(intTrajQC, 1);
			intTrajQCdet = intTrajQC - repmat(meanIntTrend, [cellNumQC 1]);
		else
			meanIntTrend = mean(intTraj, 1);
			intTrajQCdet = intTraj - repmat(meanIntTrend, [cellNum 1]);
		end
		meanInt = mean(intTrajQCdet, 2);
		intNoiseQC = intTrajQCdet - repmat(meanInt, [1 trjDuration]);
		
		% Calculate mean population intensity
		intMeanQC = mean(meanIntTrend);
		
		% Calculate zero crossing
		temp = intNoiseQC(1:end-1) .* wshift('1D', intNoiseQC(1:end-1), 1);
		intZCrossQC = sum(temp <= 0) / trjDuration / meanFrameDuration * 3600;
		
		% Calculate the ACF for intensity noise with QC
		totalAcf = zeros(1, trjDuration);
		if (useQC)
			intAcfTrajQC = zeros(cellNumQC, trjDuration);
			for k = 1:cellNumQC
				[acf, lag] = xcorr(intNoiseQC(k, :), 'biased');
				acf(lag < 0) = [];
				totalAcf = totalAcf + acf;
				intAcfTrajQC(k, :) = acf;
			end
			totalAcf = totalAcf / cellNumQC;
		else
			intAcfTrajQC = zeros(cellNum, trjDuration);
			for k = 1:cellNum
				[acf, lag] = xcorr(intNoiseQC(k, :), 'biased');
				acf(lag < 0) = [];
				totalAcf = totalAcf + acf;
				intAcfTrajQC(k, :) = acf;
			end
			totalAcf = totalAcf / cellNum;
		end
		
		% Variance (sigma2) is the zero-lag value of ACF
		intSigma2QC = totalAcf(1);
		
		% Cross correlation distance. Use squareform() to convert into
		% square matrix
		intNoiseXCorrQC = pdist(intNoiseQC, @xcorrDist);
		
		% Interpolate ACF linearly and find half maximum value
		try
			intTauHalfQC = interp1(totalAcf, 1:trjDuration, 0.5 * max(totalAcf));
			intTauHalfQC = intTauHalfQC * meanFrameDuration / 3600;
			intCV2QC = intSigma2QC / intMeanQC ^ 2;
		catch
			intTauHalfQC = -1;
			intCV2QC = -1;
		end
		
		% Check if this well has been processed for noise. If not,
		% create a new field for the well.
		if (~isfield(intNoiseInfo, wellName))
			intNoiseInfo.(wellName) = struct();
		end
		
		% Create structure to hold information for current well
		if (useQC)
			intNoiseInfo.(wellName).intMeanQC = intMeanQC;
			intNoiseInfo.(wellName).intSigma2QC = intSigma2QC;
			intNoiseInfo.(wellName).intTauHalfQC = intTauHalfQC;
			intNoiseInfo.(wellName).intCV2QC = intCV2QC;
			intNoiseInfo.(wellName).intZCrossQC = intZCrossQC;
			intNoiseInfo.(wellName).nCellQC = cellNumQC;
			intNoiseInfo.(wellName).genTrendQC = meanIntTrend;
			intNoiseInfo.(wellName).intNoiseXCorrQC = intNoiseXCorrQC;
		else
			intNoiseInfo.(wellName).intMean = intMeanQC;
			intNoiseInfo.(wellName).intSigma2 = intSigma2QC;
			intNoiseInfo.(wellName).intTauHalf = intTauHalfQC;
			intNoiseInfo.(wellName).intCV2 = intCV2QC;
			intNoiseInfo.(wellName).intZCross = intZCrossQC;
			intNoiseInfo.(wellName).nCell = cellNum;
			intNoiseInfo.(wellName).genTrend = meanIntTrend;
			intNoiseInfo.(wellName).intNoiseXCorr = intNoiseXCorrQC;
		end
		
		% Save the data to .mat file
		if (useQC)
			intNoiseTrajQC = intNoiseQC;
			save(filename, 'intMeanQC', 'intSigma2QC', 'intTauHalfQC', ...
			'intCV2QC', 'intAcfTrajQC', 'intNoiseTrajQC', 'intNoiseXCorrQC', ...
			'-append');
		else
			intMean = intMeanQC;
			intSigma2 = intSigma2QC;
			intTauHalf = intTauHalfQC;
			intCV2 = intCV2QC;
			intAcfTraj = intAcfTrajQC;
			intNoiseTraj = intNoiseQC;
			intNoiseXCorr = intNoiseXCorrQC;
			save(filename, 'intMean', 'intSigma2', 'intTauHalf', ...
			'intCV2', 'intAcfTraj', 'intNoiseTraj', 'intNoiseXCorr', '-append');
		end
		fprintf('Done!\n');
	end
	save('expInfo.mat', 'intNoiseInfo', '-append');
end