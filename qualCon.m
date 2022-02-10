% qualCon.m - Quality control
% Remove saturated trajectories
% Remove fast-changing trajectories

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

% Flag for completion of at least one round of QC
QCdone = true(1);

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

% while (1)
% 	ifRmLow = input('Remove trajectories according to max-to-population-min ratio? (Y/N): ', 's');
% 	% Detect if the length of the answer is 1, and if the answer is Y(y) or
% 	% N(n). If it is not, keep asking
% 	if (length(ifRmLow) == 1)
% 		if (ifRmLow(1) == 'y' || ifRmLow(1) == 'Y')
% 			ifRmLow = true(1);
% 			% Ask the user for the threshold of removal for trajectories
% 			% with fast changes
% 			while (1)
% 				lowThresh = input('Minimum max-to-population-min ratio for a trajectory to be included: ');
% 				% Detect if answer is a number, is at least one, and is an integer
% 				% If it is not, keep asking
% 				if (isnumeric(lowThresh) && lowThresh >= 1)
% 					break;
% 				end
% 			end
% 			break;
% 		elseif (ifRmLow(1) == 'n' || ifRmLow(1) == 'N')
% 			ifRmLow = false(1);
%			lowThresh = -1;
% 			break;
% 		end
% 	end
% end
lowThresh = -1;
ifRmLow = false(1);

if (~ifRmSat && ~ifRmFast && ~ifAcfCtrl && ~ifRmLow)
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
	
	if (ifRmLow)
		popMin = min(min(intTraj));
		for k = 1:cellNum
			if (max(intTraj(k, :)) / popMin < lowThresh)
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

fprintf('Completed QC of %d .mat file.\n', nWellProc);
save('expInfo.mat', 'ifRmSat', 'nRmSat', 'ifRmFast', 'nRmFast', ...
	'QCdone', 'ifAcfCtrl', 'acfThresh', 'ifRmLow', 'lowThresh', '-append');