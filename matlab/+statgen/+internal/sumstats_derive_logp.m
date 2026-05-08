function out = sumstats_derive_logp(pvec)
% Derive -log10(p) using the statgen sumstats zero-p convention.
    out = nan(size(pvec));
    finite_mask = isfinite(pvec);
    in_range = finite_mask & pvec >= 0 & pvec <= 1;
    zero_mask = in_range & pvec == 0;
    pos_mask = in_range & pvec > 0;
    out(zero_mask) = inf;
    out(pos_mask) = -log10(pvec(pos_mask));
end
