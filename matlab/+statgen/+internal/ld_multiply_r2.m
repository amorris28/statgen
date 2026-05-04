function out = ld_multiply_r2(ld_r2, M)
% Multiply a shard-local sparse LD-r matrix squared by an aligned payload.
    if size(M, 1) ~= size(ld_r2, 1)
        error('statgen:ld', 'ld_multiply_r2: matrix row count must match LD shard size');
    end
    out = ld_r2 * M;
end
