clear all

Nrand = 100;
times = zeros(Nrand,1);

p = gcp('nocreate');
for n_rand=1:Nrand
    spmd
        L_local = 1;
        tic
        L = gop(@plus, L_local);
        mytime = toc;
    end

    % save average from workers
    times(n_rand) = sum([mytime{:}])/p.NumWorkers;
end

% print average from trials
mean(times)

