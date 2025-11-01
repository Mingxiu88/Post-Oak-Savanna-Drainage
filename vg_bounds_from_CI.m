function [VWC_mean_tt, VWC_lo_tt, VWC_hi_tt] = vg_bounds_from_CI(MPtt, P)
    t  = MPtt.Time; vn = MPtt.Properties.VariableNames;
    nT = height(MPtt); nZ = width(MPtt);
    Vmean = nan(nT,nZ); Vlo = Vmean; Vhi = Vmean;

    for z = 1:nZ
        h = abs(MPtt{:,z});                 % suction in cm; NaNs propagate
        % --- Mean curve
        Vmean(:,z) = vg_theta(h, P.mean.thr(z), P.mean.ths(z), P.mean.alp(z), P.mean.n(z));

        % --- 95% CI envelope from 16 endpoint combinations
        Theta_all = nan(16, nT);
        idx = 0;
        for thr = [P.lo.thr(z),  P.hi.thr(z)]
          for ths = [P.lo.ths(z), P.hi.ths(z)]
            for alp = [P.lo.alp(z), P.hi.alp(z)]
              for nn = [P.lo.n(z),   P.hi.n(z)]
                  idx = idx + 1;
                  % ensure physical consistency
                  ths_use = max(ths, thr + 1e-4);
                  nn_use  = max(nn, 1.01);
                  alp_use = max(alp, 1e-6);
                  Theta_all(idx,:) = vg_theta(h, thr, ths_use, alp_use, nn_use);
              end
            end
          end
        end
        Vlo(:,z) = min(Theta_all, [], 1, 'omitnan').';
        Vhi(:,z) = max(Theta_all, [], 1, 'omitnan').';
    end

    VWC_mean_tt = array2timetable(Vmean,'RowTimes',t,'VariableNames',vn);
    VWC_lo_tt   = array2timetable(Vlo,  'RowTimes',t,'VariableNames',vn);
    VWC_hi_tt   = array2timetable(Vhi,  'RowTimes',t,'VariableNames',vn);
end