function theta = vg_theta(h_cm, thr, ths, alpha, nvg)
    m = 1 - 1./nvg;
    Se = (1 + (alpha .* h_cm).^nvg) .^ (-m);
    theta = thr + (ths - thr) .* Se;
    % numerical guard
    theta = min(max(theta, thr), ths);
end