function out = fast_prune(logpvec, ld_panel, r2_threshold, chrX_sex)
%FAST_PRUNE Greedy significance-based LD pruning.
%
%   out = statgen.fast_prune(logpvec, ld_panel)
%   out = statgen.fast_prune(logpvec, ld_panel, r2_threshold)
%   out = statgen.fast_prune(logpvec, ld_panel, r2_threshold, chrX_sex)
%
% logpvec is a num_snp vector aligned to ld_panel. The output keeps retained
% scores and sets pruned neighbors to NaN. r2_threshold defaults to 0.2.
% chrX_sex selects the chrX LD shard when chrX is loaded.
%
% See also statgen.LDPanel, statgen.LDPanel.multiply_r2.
    if nargin < 3 || isempty(r2_threshold)
        r2_threshold = 0.2;
    end
    if nargin < 4
        chrX_sex = [];
    end
    if ~isscalar(r2_threshold) || ~isfinite(r2_threshold) || r2_threshold < 0
        error('statgen:ld', 'fast_prune: r2_threshold must be a finite non-negative number');
    end

    original_size = size(logpvec);
    row_vector = isrow(logpvec) && numel(logpvec) == ld_panel.num_snp;
    if ~isvector(logpvec)
        error('statgen:ld', 'fast_prune: logpvec must be a vector');
    end
    values = double(logpvec(:));
    if numel(values) ~= ld_panel.num_snp
        error('statgen:ld', ...
            'fast_prune: vector length mismatch: expected %d, got %d', ...
            ld_panel.num_snp, numel(values));
    end

    selected = statgen.internal.ld_selected_shards(ld_panel, chrX_sex);
    out = values;
    for i = 1:numel(selected)
        offset = ld_panel.shard_offsets(i);
        rows = (offset.start0 + 1):offset.stop0;
        out(rows) = prune_shard_(out(rows), selected{i}.ld_r2, r2_threshold);
    end

    if row_vector
        out = reshape(out, original_size);
    end
end

function out = prune_shard_(values, ld_r2, threshold)
    out = values(:);
    finite_idx = find(isfinite(out));
    if isempty(finite_idx)
        return
    end

    keys = [-abs(out(finite_idx)), finite_idx(:)];
    [~, order_idx] = sortrows(keys, [1 2]);
    order = finite_idx(order_idx);

    pruned = false(numel(out), 1);
    retained = false(numel(out), 1);
    for ii = 1:numel(order)
        idx = order(ii);
        if pruned(idx)
            continue
        end
        retained(idx) = true;
        [rows, ~, vals] = find(ld_r2(:, idx));
        if isempty(rows)
            continue
        end
        neighbors = rows(vals >= threshold);
        neighbors = neighbors(neighbors ~= idx);
        neighbors = neighbors(~retained(neighbors));
        if ~isempty(neighbors)
            pruned(neighbors) = true;
            out(neighbors) = NaN;
        end
    end
end
