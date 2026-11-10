function extract_svc_for_analysis(SPMmat, mask_file, contrasts, threshold, extent)
% EXTRACT_SVC_FOR_ANALYSIS  Small-volume correction (RFT) within a mask.
%
% extract_svc_for_analysis(SPMmat, mask_file, contrasts, threshold, extent)
%
% Reproduces SPM's "small volume correction" without the GUI: resel counts
% are computed inside the mask (spm_resels_vol) and the FWE peak threshold
% from random field theory (spm_uc_RF).
%
% Per contrast index in `contrasts`:
%   1. threshold at p < .001 uncorrected (spm_getSPM)
%   2. keep suprathreshold voxels inside mask_file
%   3. compute the SVC FWE threshold at `threshold` (e.g. .05)
%   4. if any voxel survives, write to the SPM.mat folder:
%        SVC_results_<con>.mat / .csv   (coordinates and Z of survivors)
%        spmT_NNNN_SVC_significant_voxels_<con>_FWE_p05.nii  (binary mask)

spmDir = fileparts(SPMmat);
load(SPMmat, 'SPM');
cd(spmDir);

% Read the mask once
Vm        = spm_vol(mask_file);
mask_data = spm_read_vols(Vm) > 0;

for iCon = contrasts

    if iCon > numel(SPM.xCon)
        fprintf('  Contrast %d does not exist in SPM.mat (only %d defined).\n', ...
            iCon, numel(SPM.xCon));
        continue;
    end

    conName = SPM.xCon(iCon).name;
    fprintf('\n  Processing contrast %d: %s\n', iCon, conName);

    % ---- Threshold at p < 0.001 uncorrected ---------------------------
    xSPM = struct();
    xSPM.swd       = spmDir;
    xSPM.title     = sprintf('SVC_%s', conName);
    xSPM.Ic        = iCon;
    xSPM.n         = 0;
    xSPM.Im        = [];
    xSPM.pm        = [];
    xSPM.Ex        = [];
    xSPM.u         = 0.001;
    xSPM.k         = extent;
    xSPM.thresDesc = 'none';

    try
        [~, xSPM] = spm_getSPM(xSPM);

        if isempty(xSPM.XYZ)
            fprintf('    No voxels at p < 0.001 uncorrected.\n');
            continue;
        end
        fprintf('    Found %d voxels at p < 0.001 uncorrected.\n', size(xSPM.XYZ, 2));

        % ---- Restrict to voxels inside the mask -----------------------
        in_mask = false(1, size(xSPM.XYZ, 2));
        for i = 1:size(xSPM.XYZ, 2)
            vox = round(Vm.mat \ [xSPM.XYZmm(:,i); 1]);
            vox = vox(1:3);
            if all(vox >= 1) && all(vox(:)' <= Vm.dim) ...
                    && mask_data(vox(1), vox(2), vox(3))
                in_mask(i) = true;
            end
        end

        if ~any(in_mask)
            fprintf('    No voxels inside the mask.\n');
            continue;
        end
        fprintf('    %d voxels inside the mask.\n', sum(in_mask));

        % ---- Compute resel counts within the mask ---------------------
        % spm_resels_vol(V, FWHM) returns 4-element resel count R
        % where R(1) = volume in resels, etc. FWHM comes from SPM.xVol.
        R_svc = spm_resels_vol(Vm, SPM.xVol.FWHM)';

        % Degrees of freedom for this contrast
        df_con = [SPM.xCon(iCon).eidf, SPM.xX.erdf];

        % STAT type ('T' or 'F')
        STAT = SPM.xCon(iCon).STAT;

        % ---- RFT FWE peak threshold within the small volume -----------
        u_FWE = spm_uc_RF(threshold, df_con, STAT, R_svc, 1);
        fprintf('    RFT SVC FWE threshold: %s = %.3f (df = [%g %g], resels = %.2f)\n', ...
            STAT, u_FWE, df_con(1), df_con(2), R_svc(end));

        % ---- Find in-mask voxels exceeding the threshold --------------
        survives = in_mask & (xSPM.Z >= u_FWE);

        if ~any(survives)
            fprintf('    No in-mask voxels survive RFT SVC at p < %.3f FWE.\n', threshold);
            continue;
        end
        fprintf('    %d in-mask voxels survive RFT SVC at p < %.3f FWE.\n', ...
            sum(survives), threshold);

        % ---- Save: .mat, .csv, and NIfTI mask -------------------------
        results              = struct();
        results.contrastName = conName;
        results.mask         = mask_file;
        results.threshold    = threshold;
        results.u_FWE        = u_FWE;
        results.R_svc        = R_svc;
        results.XYZ          = xSPM.XYZ(:, survives);
        results.XYZmm        = xSPM.XYZmm(:, survives);
        results.Z            = xSPM.Z(survives);

        save(fullfile(spmDir, sprintf('SVC_results_%s.mat', conName)), 'results');

        T = table(results.XYZmm(1,:)', results.XYZmm(2,:)', ...
                  results.XYZmm(3,:)', results.Z', ...
                  'VariableNames', {'X_mm','Y_mm','Z_mm','Z_score'});
        writetable(T, fullfile(spmDir, sprintf('SVC_results_%s.csv', conName)));

        % NIfTI mask of surviving voxels, sized to the spmT volume
        Vref    = spm_vol(fullfile(spmDir, sprintf('spmT_%04d.nii', iCon)));
        maskout = zeros(Vref.dim);
        linidx  = sub2ind(Vref.dim, results.XYZ(1,:)', ...
                                    results.XYZ(2,:)', ...
                                    results.XYZ(3,:)');
        maskout(linidx) = 1;

        Vout         = Vref;
        Vout.fname   = fullfile(spmDir, sprintf( ...
            'spmT_%04d_SVC_significant_voxels_%s_FWE_p05.nii', iCon, conName));
        Vout.dt      = [spm_type('uint8') 0];
        Vout.descrip = sprintf('SVC RFT FWE p<%.3f, %s', threshold, conName);
        Vout = rmfield(Vout, 'pinfo');
        spm_write_vol(Vout, maskout);

        fprintf('    Saved .mat, .csv, and NIfTI mask.\n');

    catch ME
        fprintf('    SVC failed for contrast %d (%s): %s\n', iCon, conName, ME.message);
    end
end

end