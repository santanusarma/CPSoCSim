function list_temperatures(T, N_core_x, N_blocks_per_core)

fprintf('Core\tBlock\tTemperature (C)\n\n');
i = 1;
j = 1;
k = 1;
for q = 1:length(T)
    fprintf('(%d,%d)\t%d\t%.1f\n', i, j, k, T(q));
    k = k + 1;
    if k > N_blocks_per_core
        k = 1;
        j = j + 1;
    end
    if j > N_core_x
        j = 1;
        i = i + 1;
    end
end
