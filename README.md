# Noise Analysis Breakdown

In this documentation, a step-by-step description of the code used in the noise anlaysis of the HIV Gene Expression assay.

Lu, Y., Bohn-Wippert, K., Pazerunas, P. J., Moy, J. M., Singh, H., & Dar, R. D. (2021). Screening for gene expression fluctuations reveals latency-promoting agents of HIV. Proceedings of the National Academy of Sciences of the United States of America, 118(11), e2012191118. https://doi.org/10.1073/pnas.2012191118



***
### 1) Read .mat files from xlsx2mat.m and Plot Trajectories for Each Well
***

*Read .mat files & plot trajectories for each well's raw data. Optionally, plot equivalent diameter and save plots as .png files.*

<details>
	<summary><b>1.First we check to see if file has already been read and plotted.</b></summary>

```python
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
```
  </details>
  
<details>
	<summary><b>2. Ask user if equivalent diameter should be plotted.</b></summary>

```python
while (1)
	ifEqDiam = input('Do you want to plot equivalent diameter? (Y/N): ', 's');
	% Detect if the length of the answer is 1, and if the answer is Y(y) or
	% N(n). If it is not, keep asking
	if (length(ifEqDiam) == 1)
		if (ifEqDiam(1) == 'y' || ifEqDiam(1) == 'Y')
			ifEqDiam = true(1);
			break;
		elseif (ifEqDiam(1) == 'n' || ifEqDiam(1) == 'N')
			ifEqDiam = false(1);
			break;
		end
	end
end
```
  </details>
  
<details>
	<summary><b>3. Ask user to save plot as .png</b></summary>

```python
while (1)
	ifSave = input('Do you want to save the plots? (Y/N): ', 's');
	% Detect if the length of the answer is 1, and if the answer is Y(y) or
	% N(n). If it is not, keep asking
	if (length(ifSave) == 1)
		if (ifSave(1) == 'y' || ifSave(1) == 'Y')
			ifSave = true(1);
			break;
		elseif (ifSave(1) == 'n' || ifSave(1) == 'N')
			ifSave = false(1);
			break;
		end
	end
end
```
  </details>
  
<details>
	<summary><b>4. Ask user to include general trends.</b></summary>
  
```python
while (1)
	ifGenTrend = input('Do you want to include general trends? (Y/N): ', 's');
	% Detect if the length of the answer is 1, and if the answer is Y(y) or
	% N(n). If it is not, keep asking
	if (length(ifGenTrend) == 1)
		if (ifGenTrend(1) == 'y' || ifGenTrend(1) == 'Y')
			ifGenTrend = true(1);
			break;
		elseif (ifGenTrend(1) == 'n' || ifGenTrend(1) == 'N')
			ifGenTrend = false(1);
			break;
		end
	end
end
```
  </details>

<details>
	<summary><b>5. User can re-plot data after running Quality Control (Step 4) function.</b></summary>
  
```python
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
```
  </details>

<details>
	<summary><b>6. Plot Data!</b></summary>
  
```python
for i = 1:nWellProc
	% Write filename
	filename = sprintf('Well%s%s.mat', char(rowList(i)), char(colList(i)));

	% Load data
	load(filename);
	fprintf('Reached well #%d/%d! Plotting... ', i, nWellProc);
	

	% Create a new figure and plot
	f = figure;
	hold on;
	grid on;
	if (cellNum == 0)
		fprintf('No cells detected!\n');
		if (ifSave)
			if (exist('Trajectory Plots', 'dir') ~= 7)
				mkdir('Trajectory Plots');
			end
			if (exist('tag', 'var') == 1 && ~isempty(tag))
				if (exist(sprintf('Trajectory Plots\\%s', tag), 'dir') ~= 7)
					mkdir(sprintf('Trajectory Plots\\%s', tag));
				end
				saveas(f, sprintf('Trajectory Plots\\%s\\Intensity Well%s%s.png', ...
					tag, char(rowList(i)), char(colList(i))));
			else
				saveas(f, sprintf('Trajectory Plots\\Intensity Well%s%s.png', ...
					char(rowList(i)), char(colList(i))));
			end
			close(f);
			drawnow;
		end
		continue;
	end
	
	if (exist('tag', 'var') == 1 && ~isempty(tag))
		tempString = sprintf('%s%s: %s', char(rowList(i)), char(colList(i)), tag);
	else
		tempString = sprintf('%s%s', char(rowList(i)), char(colList(i)));
	end
	
	tPlot = ((trjStart - 1):(trjDuration + trjStart - 2)) / 3600 * meanFrameDuration;
	
	if (ifQC)
		dIntTraj = intTrajQC - wshift('2d', intTrajQC, [0 -1]);
		for k = 1:cellNumQC
			plot(tPlot, intTrajQC(k, :));
		end
		for k = 2:trjDuration
			sigmaDInt(k) = std(dIntTraj(:, k));
			muDInt(k) = mean(dIntTraj(:, k));
		end
		if (ifGenTrend)
			genTrend = mean(intTrajQC, 1);
			plot(tPlot, genTrend, 'k', 'LineWidth', 2);
		end
		yyaxis right;
		hold on;
		if (cellNum >= 200)
			thresh = 200 / cellNumQC;
		else
			thresh = 2 - 1 / 200 * cellNumQC;
		end
		plot(tPlot, sigmaDInt./abs(muDInt)<=thresh, 'b', 'LineWidth', 2)
		ylim([0 10])
		ax = gca;
		ax.YColor = [0 0 1];
		
		title(sprintf('%d intensity trajectories from well %s', ...
			cellNumQC, tempString));
	else
		dIntTraj = intTraj - wshift('2d', intTraj, [0 -1]);
		for k = 1:cellNum
			plot(tPlot, intTraj(k, :));
		end
		for k = 2:trjDuration
			sigmaDInt(k) = std(dIntTraj(:, k));
			muDInt(k) = mean(dIntTraj(:, k));
		end
		if (ifGenTrend)
			genTrend = mean(intTraj, 1);
			plot(tPlot, genTrend, 'k', 'LineWidth', 2);
		end
		yyaxis right;
		hold on;
		if (cellNum >= 200)
			thresh = 200 / cellNum;
		else
			thresh = 2 - 1 / 200 * cellNum;
		end
		plot(tPlot, sigmaDInt./abs(muDInt)<=thresh, 'b', 'LineWidth', 2)
		ylim([0 10])
		ax = gca;
		ax.YColor = [0 0 1];
		
		title(sprintf('%d intensity trajectories from well %s', ...
			cellNum, tempString));
	end
	xlim([(trjStart - 1) (trjStart + trjDuration - 2)] * meanFrameDuration / 3600);
	xlabel('Time (h)');
	yyaxis left;
	ylabel('Intensity (a.u.)');

	% Save the figures if the user required so
	if (ifSave)
		if (exist('Trajectory Plots', 'dir') ~= 7)
			mkdir('Trajectory Plots');
		end
		if (exist('tag', 'var') == 1 && ~isempty(tag))
			if (exist(sprintf('Trajectory Plots\\%s', tag), 'dir') ~= 7)
				mkdir(sprintf('Trajectory Plots\\%s', tag));
			end
			saveas(f, sprintf('Trajectory Plots\\%s\\Intensity Well%s%s.png', ...
				tag, char(rowList(i)), char(colList(i))));
		else
			saveas(f, sprintf('Trajectory Plots\\Intensity Well%s%s.png', ...
				char(rowList(i)), char(colList(i))));
		end
		close(f);
		drawnow;
	end

	% Plot equivalent diameter if the user required so
	if (ifEqDiam)
		f = figure;
		hold on;
		grid on;
		for k = 1:cellNum
			plot(timeTraj(k, :) / 3600, eqDiamTraj(k, :));
		end
		xlim([(trjStart - 1) (trjStart + trjDuration - 1)] * meanFrameDuration / 3600);
		xlabel('Time (h)');
		ylabel('Equivalent Diameter (\mum)');
		title(sprintf('%d equivalent diameter trajectories from well %s', ...
			cellNum, tempString));
		
		if (ifSave)
			if (exist('Trajectory Plots', 'dir') ~= 7)
				mkdir('Trajectory Plots');
			end
			if (exist('tag', 'var') == 1 && ~isempty(tag))
				if (exist(sprintf('Trajectory Plots\\%s', tag), 'dir') ~= 7)
					mkdir(sprintf('Trajectory Plots\\%s', tag));
				end
				saveas(f, sprintf('Trajectory Plots\\%s\\EqDiam Well%s%s.png', ...
					tag, char(rowList(i)), char(colList(i))));
			else
				saveas(f, sprintf('Trajectory Plots\\EqDiam Well%s%s.png', ...
					char(rowList(i)), char(colList(i))));
			end
			close(f);
			drawnow;
		end
	end

	% Report progress
	fprintf('Completed!\n');
end

fprintf('Completed plotting of %d .mat file!\n', nWellProc);
```
   </details>
  

### 2) Noise Analysis for Raw Intensity Trajectory Data
***
* The noise analysis was first run over the raw data (before auto-correlation in Quality Control). 
* Can only be used over isoclonal data
* Requires Step 1 

<details>
	<summary><b>1. Ensure .mat file has been read and plotted.</b></summary>
  
```python
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
```
  </details>
  
<details>
  <summary>2. Check to see if noise analysis has already been performed. If so, we want to add to the previous file. This is useful when we re-run the noise analysis over the auto-correlated data.</summary>
  
```python
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
```
  </details>

<details>
  <summary>3. Check to see if auto-correlation (Quality Control) has been perfromed over data. Asks user's input on whether or not to use QC data. Since we have not yet performed Quality Control over the data, we will proceed with the NonQC process for now.</summary>
  
```python
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
```
  </details>

<details>
  <summary>4. Set QC condition based on previous step.</summary>
  
```python
if (useQC)
	intNoiseInfo.existQC = true;
else
	intNoiseInfo.existNonQC = true;
end
```
  </details>

<details>
  <summary>5. Perform noise analysis! </summary>
	
```python
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
```					    
  </details>

### 3) Quality Control
*****
* Includes an auto-correlation fucntion
* Removes saturated and fast-changing trajectories
* Requires Step 1 & 3

<details>
  <summary>1. Check to see if .mat file has been read and raw data trajectories have been plotted (Step 1).</summary>
  
```python
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

% Load data
load('expInfo.mat');
```
  </details>
  
<details>
  <summary>2. Label data as having been auto-correlated.</summary>
  
```python
% Flag for completion of at least one round of QC
QCdone = true(1);
```
  </details>

<details>
  <summary>3. Ask user if they would like to remove saturated trajectories.</summary>
  
```python
% Ask the user whether to remove trajectories with long saturations
while (1)
	ifRmSat = input('Remove trajectories with long saturations? (Y/N): ', 's');
	% Detect if the length of the answer is 1, and if the answer is Y(y) or
	% N(n). If it is not, keep asking
	if (length(ifRmSat) == 1)
		if (ifRmSat(1) == 'y' || ifRmSat(1) == 'Y')
			ifRmSat = true(1);
			% Ask the user for the threshold of removal for trajectories
			% with saturations
			while (1)
				nRmSat = input('Number of frames at saturation for a trajectory to be excluded: ');
				% Detect if answer is a number, is at least one, and is an integer
				% If it is not, keep asking
				if (isnumeric(nRmSat) && nRmSat >= 1 && floor(nRmSat) == nRmSat)
					break;
				end
			end
			break;
		elseif (ifRmSat(1) == 'n' || ifRmSat(1) == 'N')
			ifRmSat = false(1);
			nRmSat = -1;
			break;
		end
	end
end
```
  </details>

<details>
  <summary>4. Ask user if they would like to remove fast-changing trajectories.</summary>
  
```python
% Ask user whether to remove fast changing trajectories
while (1)
	ifRmFast = input('Remove trajectories that changes too fast? (Y/N): ', 's');
	% Detect if the length of the answer is 1, and if the answer is Y(y) or
	% N(n). If it is not, keep asking
	if (length(ifRmFast) == 1)
		if (ifRmFast(1) == 'y' || ifRmFast(1) == 'Y')
			ifRmFast = true(1);
			% Ask the user for the threshold of removal for trajectories
			% with fast changes
			while (1)
				nRmFast = input('Maximum percent change of fluorescence between adjacent frames for \na trajectory to be included: ');
				% Detect if answer is a number, is at least one, and is an integer
				% If it is not, keep asking
				if (isnumeric(nRmFast) && nRmFast >= 1 && floor(nRmFast) == nRmFast)
					break;
				end
			end
			break;
		elseif (ifRmFast(1) == 'n' || ifRmFast(1) == 'N')
			ifRmFast = false(1);
			nRmFast = -1;
			break;
		end
	end
end
```
  </details>

<details>
  <summary>5. Ask user if they would like to remove trajectories with siginificant white noise.</summary>
  
```python
% Ask user whether to remove trajectories with large white noise according
% to ACF
if (trjDuration >= 2)
	while (1)
		ifAcfCtrl = input('Remove trajectories with significant white noise? (Y/N): ', 's');
		% Detect if the length of the answer is 1, and if the answer is Y(y) or
		% N(n). If it is not, keep asking
		if (length(ifAcfCtrl) == 1)
			if (ifAcfCtrl(1) == 'y' || ifAcfCtrl(1) == 'Y')
				ifAcfCtrl = true(1);
				% Ask the user for the threshold of removal for trajectories
				% with fast changes
				while (1)
					acfThresh = input('Maximum percentage drop of ACF from first frame to the second for \na trajectory to be included: ');
					% Detect if answer is a number, is at least one, and is an integer
					% If it is not, keep asking
					if (isnumeric(acfThresh) && acfThresh <= 100 && acfThresh >= 0)
						break;
					end
				end
				break;
			elseif (ifAcfCtrl(1) == 'n' || ifAcfCtrl(1) == 'N')
				ifAcfCtrl = false(1);
				acfThresh = -1;
				break;
			end
		end
	end
end
```
  </details>

<details>
  <summary>6. Auto-correlate data</summary>
  
```python
lowThresh = -1;
ifRmLow = false(1);

if (~ifRmSat && ~ifRmFast && ~ifAcfCtrl)
	return;
end

for i = 1:nWellProc
	fprintf('Reached well #%d/%d. Processing... ', i, nWellProc);
	% Write filename
	filename = sprintf('Well%s%s.mat', char(rowList(i)), char(colList(i)));
	
	% Load data
	load(filename);
	
	if (cellNum == 0)
		fprintf('No cells detected!\n');
		continue;
	end

	% Initialize array for removal
	rmIndex = false(1, cellNum);

	% Mark saturated trajectories
	if (ifRmSat)
		% Step through all trajectories
		for k = 1:cellNum
			% Variable for storing number of consecutive saturation
			% frames, and the max value of it in a single trajectory
			curSatLength = 0;
			maxSatLength = 0;
			% Step through every frame of a trajectory
			for l = 1:trjDuration
				% If a frame crosses the threshold or not. 
				% 65535 is the maximum intensity reported by the camera
				% The tolerance is 95% of the maximum intensity
				% If the frame exceeds 95% of maximum, it is considered
				% saturated.
				if (intTraj(k, l) >= 0.95 * 65535)
					curSatLength = curSatLength + 1;
				else
					% If a frame does not exceed 95% of maximum, check
					% if we were tracking a saturated segment. If so,
					% compare and store the larger value between the
					% recorded maximum length, and the current length,
					% then reset the current length to 0.
					if (curSatLength ~= 0)
						maxSatLength = max(maxSatLength, curSatLength);
						curSatLength = 0;
					end
				end
			end
			% If the maximum saturation length exceeds user-assigned
			% threshold, mark it for exclusion
			if (maxSatLength >= nRmSat)
				rmIndex(k) = true(1);
			end
		end
	end

	% Mark fast-changing trajectories
	if (ifRmFast)
		% Step through all trajectories
		for k = 1:cellNum
			dInt = abs(intTraj(k, :) - wshift(1, intTraj(k, :), -1));
			dInt(1) = [];
			dInt = dInt ./ intTraj(k, 1:(end - 1));
			if (sum(dInt >= nRmFast/100) > 0)
				rmIndex(k) = true(1);
			end
		end
	end
	
	if (ifAcfCtrl)
		for k = 1:cellNum
			if (intAcfTraj(k, 2) / intAcfTraj(k, 1) <= 1 - acfThresh / 100)
				rmIndex(k) = true;
			end
		end
	end
	
	% Remove trajectories
	eqDiamTrajQC = eqDiamTraj(rmIndex == 0, :);
	intTrajQC = intTraj(rmIndex == 0, :);
	timeTrajQC = timeTraj(rmIndex == 0, :);
	timePtTrajQC = timePtTraj(rmIndex == 0, :);
	cellNumQC = cellNum - sum(rmIndex);

	% Append QC data to .mat file
	save(filename, 'eqDiamTrajQC', 'intTrajQC', 'timeTrajQC', ...
		'timePtTrajQC', 'cellNumQC', '-append');

	% Report progress
	fprintf('Completed!\n%d trajectories removed from %s.\n\n', ...
		sum(rmIndex), filename);
end
```
  </details>
  
### 4) Repeat Noise Analysis for QC Intensity Trajectory Data
****
* Repeated after auto-correlation in Quality Control step
* Same process as Step 2
* QC data is considered the "final data." Thus, this process is the main step in the noise analysis.


### 5) Add Tags to WellName.mat files
*****
* Label wells as positive/negative control and data points

<details>
  <summary>1. Check to see if .mat files have been read (Step 1).</summary>
  
```python
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
```
  </details>

<details>
  <summary>2. User can tag wells using an Excel sheet containing information on how to tag the wells. If not, user can manually add tagging data later on.</summary>
	
```python
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
```
  </details>

<details>
	<summary>3. Create a list with tags. If Excel data not used, prompt user to manually add tags for each file.</summary>
	
```python
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
```
</details>
	
<details>
	<summary>4. Iterate over tag list and count number of unique tags.</summary>
	
```python
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
```
</details>

### 6) Exclude Wells
****
* User can decide to remove any wells from analysis. This applies to empty wells, data with significant systematic shift, or specific wells that the user can manually input. 

<details>
  <summary>1. Check to see if .mat files have been read (Step 1).</summary>
  
```python
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
```
  </details>

<details>
  <summary>2. Ask user if they would like to use information from an Excel sheet to remove sepcific wells. </summary>
  
```python
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
```
  </details>

<details>
  <summary>3. User can remove wells with data that varies significantly from expected values. </summary>
  
```python
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
```
  </details>

<details>
  <summary>4. User can remove wells that contain no data.</summary>
  
```python
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
```
  </details>

<details>
  <summary>5. If user opts not to use an Excel sheet to remove wells, they can manually input tags to remove asscociated wells. </summary>
  
```python
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
```
  </details>

<details>
  <summary>6. Removal of wells with significant systematic shift, if user opted in. </summary>
  
```python
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
```
  </details>

<details>
  <summary>7. After removal of wells, update list of wells' tags.</summary>
  
```python
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
```
  </details>

### 7) Rescan .mat files
*****
* User can decide to rescan if any wells were excluded

<details>
  <summary>1. Check to see if .mat files have been read (Step 1).</summary>
  
```python
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
```
  </details>

<details>
  <summary>2. Iterate over all wells and see if well has been saved. Save previously excluded wells and prompt user to re-run tag wells code to tag newly saved wells.</summary>
  
```python
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
```
  </details>
