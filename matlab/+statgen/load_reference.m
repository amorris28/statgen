function panel = load_reference(path, shards)
% Load a reference panel from a .bim file or sharded .bim template.
%
%   panel = statgen.load_reference(path)
%   panel = statgen.load_reference(path, shards)
%
% If path contains '@', it is a sharded template (e.g. 'chr@.bim') and
% shard discovery uses canonical-label substitution only (1-22, X).
% Otherwise path is a single file split by chr column into canonical shards.
    if nargin < 2
        shards = [];
    end

    path = char(path);
    if ~isempty(strfind(path, '@'))
        panel = load_sharded_(path, shards);
    else
        bim = reference_bim_rows_(statgen.internal.bfile_parse_bim(path));
        panel = load_split_by_chr_(bim, shards);
    end
end

% ---------------------------------------------------------------------------

function panel = load_sharded_(path, requested_shards)
    canonical = statgen.internal.canonical_labels();
    available_labels = {};
    available_paths = {};
    for i = 1:numel(canonical)
        label = canonical{i};
        bim_path = strrep(path, '@', label);
        if exist(bim_path, 'file') == 2
            available_labels{end+1} = label; %#ok<AGROW>
            available_paths{end+1} = bim_path; %#ok<AGROW>
        end
    end

    if isempty(available_labels)
        error('statgen:io', 'No BIM shards found matching template: %s', path);
    end

    selected = statgen.internal.validate_requested_shards(requested_shards, available_labels, 'load_reference');
    shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        path_idx = find(strcmp(available_labels, label), 1, 'first');
        bim = reference_bim_rows_(statgen.internal.bfile_parse_bim(available_paths{path_idx}));
        shards{i} = statgen.ReferenceShard(label, bim.chr, bim.snp, bim.bp, bim.a1, bim.a2);
    end
    panel = statgen.ReferencePanel(shards);
end

function panel = load_split_by_chr_(bim, requested_shards)
    canonical = statgen.internal.canonical_labels();
    available = {};
    for i = 1:numel(canonical)
        c = canonical{i};
        if any(strcmp(bim.chr, c))
            available{end+1} = c; %#ok<AGROW>
        end
    end

    selected = statgen.internal.validate_requested_shards(requested_shards, available, 'load_reference');
    shards = cell(numel(selected), 1);
    for ci = 1:numel(selected)
        c = selected{ci};
        mask = strcmp(bim.chr, c);
        shards{ci} = statgen.ReferenceShard(c, ...
            bim.chr(mask), bim.snp(mask), ...
            bim.bp(mask), bim.a1(mask), bim.a2(mask));
    end

    panel = statgen.ReferencePanel(shards);
end

function out = reference_bim_rows_(bim)
    canonical = statgen.internal.canonical_labels();
    keep = ismember(bim.chr, canonical);
    out.chr = bim.chr(keep);
    out.snp = bim.snp(keep);
    out.bp = bim.bp(keep);
    out.a1 = bim.a1(keep);
    out.a2 = bim.a2(keep);
end
