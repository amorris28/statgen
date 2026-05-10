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
        function obj = AnnotationShard(label, reference_checksum, annomat, validate_binary)
            if nargin == 0, return; end
            if nargin < 4 || isempty(validate_binary), validate_binary = true; end
            obj.label = char(label);
            obj.reference_checksum = char(reference_checksum);
            obj.annomat = ensure_sparse_binary_(annomat, validate_binary);
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
            fprintf('    annomat: %d-by-%d sparse logical-equivalent, nnz=%d, density=%.4g\n', ...
                obj.num_snp, obj.num_annot, nz, nz / denom);
            fprintf('    reference_checksum: %s\n', obj.reference_checksum);
        end
    end
end

function out = ensure_sparse_binary_(x, validate_binary)
    if issparse(x)
        out = sparse(x);
    else
        out = sparse(double(x));
    end
    if ~validate_binary
        return
    end
    vals = nonzeros(out);
    if any((vals ~= 0) & (vals ~= 1))
        bad = vals(find((vals ~= 0) & (vals ~= 1), 1, 'first'));
        error('statgen:annotations', 'annomat contains non-binary value: %g', bad);
    end
    out = spones(out);
end
