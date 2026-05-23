classdef AnnotationShard
% Immutable annotation matrix for one reference shard.
    properties (SetAccess = private)
        label
        reference_checksum
        num_snp
        num_annot
        annomat
    end

    methods
        function obj = AnnotationShard(label, reference_checksum, annomat)
            if nargin == 0, return; end
            obj.label = char(label);
            obj.reference_checksum = char(reference_checksum);
            obj.annomat = ensure_sparse_numeric_(annomat);
            obj.num_snp = size(obj.annomat, 1);
            obj.num_annot = size(obj.annomat, 2);
        end

        function display(obj)
            name = inputname(1);
            if ~isempty(name)
                fprintf('%s =\n\n', name);
            end
            disp(obj);
        end

        function disp(obj)
            if numel(obj) ~= 1
                fprintf('  statgen.AnnotationShard array with size %s\n', statgen.internal.display_size_string(size(obj)));
                return;
            end
            nz = nnz(obj.annomat);
            denom = max(1, obj.num_snp * obj.num_annot);
            fprintf('  statgen.AnnotationShard object\n\n');
            fprintf('    label: %s\n', obj.label);
            fprintf('    num_snp: %d\n', obj.num_snp);
            fprintf('    num_annot: %d\n', obj.num_annot);
            fprintf('    annomat: %d-by-%d sparse numeric, nnz=%d, density=%.4g\n', ...
                obj.num_snp, obj.num_annot, nz, nz / denom);
            fprintf('    reference_checksum: %s\n', obj.reference_checksum);
        end
    end
end

function out = ensure_sparse_numeric_(x)
    if issparse(x)
        out = sparse(double(x));
    else
        out = sparse(double(x));
    end
    vals = nonzeros(out);
    if any(~isfinite(vals))
        bad = vals(find(~isfinite(vals), 1, 'first'));
        error('statgen:annotations', 'annomat contains non-finite value: %g', bad);
    end
end
