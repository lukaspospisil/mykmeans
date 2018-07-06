clear all

%poolobj = parpool(2);

tic
parfor t=1:100
    pause(0.0001)
end
toc

