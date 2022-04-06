function [patient_HC_correls,spectral_deviation_index] = compute_spectraldeviationindex(patient_data,control_data)

%PURPOSE:           Compute spectral deviation index as in: "Aberrant neurophysiological signaling underlies speech impairments in Parkinson’s disease" (Wiesman et al., 2022)
%
%REQUIRED INPUTS:   patient_data: array of data for patient participants - dimensions: locations (e.g., vertices or ROIs) x participants x frequency samples
%                   control_data: array of data for control/reference participants - dimensions: locations x participants x frequency samples
%
%OUTPUTS:           patient_HC_correls: cross-correlogram (Pearson r) of power values over all frequencies between each combination of patient and control participant - dimensions: locations x patientn x controln
%                   spectral_deviation_index: deviation indices for each patient per location - dimensions: location x patient
%
%AUTHOR:            Alex I. Wiesman, neuroSPEED lab, Montreal Neurological Institute
%VERSION HISTORY:   04/06/2022  v1: First working version of program
%
%LICENSE:           This software is distributed under the terms of the GNU General Public License as published by the Free Software Foundation. Further details on the GPLv3 license can be found at http://www.gnu.org/copyleft/gpl.html.
%                   FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE AUTHORS DO NOT MAKE ANY WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.

%check if matrices are equivalently sized (dimensions 1 and 3)
if ~isequal(size(patient_data,1),size(control_data,1)) || ~isequal(size(patient_data,3),size(control_data,3))
    error('Patient and control data matrices have unequal dimensions.')
end

%compute the similarity matrix between each patient and all controls
patient_HC_correls = nan(size(control_data,1),size(patient_data,2),size(control_data,2));
for i = 1:size(control_data,2)
    fprintf('Computing correlations for HC # %d...\n',i);
    for ii = 1:size(patient_data,2)
        HC_vec = squeeze(control_data(:,i,:));
        PD_vec = squeeze(patient_data(:,ii,:));
        for iii = 1:size(HC_vec,1)
            correl = corrcoef(HC_vec(iii,:),PD_vec(iii,:));
            patient_HC_correls(iii,ii,i) = correl(1,2);
        end
    end
end

spectral_deviation_index = 1-median(atanh(patient_HC_correls),3); 