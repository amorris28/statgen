function sumstats_warn_optional_zn_completeness(zvec, nvec, logpvec, context)
% Warn when optional z/n payloads are absent or incomplete for present sumstats.
    is_present = ~isnan(logpvec);
    warn_one_optional_(zvec, 'zvec', is_present, context);
    warn_one_optional_(nvec, 'nvec', is_present, context);
end

function warn_one_optional_(vec, name, is_present, context)
    if isempty(vec)
        warning('statgen:sumstats', '%s: %s is absent', context, name);
    elseif any(isnan(vec(is_present)))
        warning('statgen:sumstats', '%s: %s has missing values among present sumstats variants', context, name);
    end
end
