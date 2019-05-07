function [phat,cint] = binofit_wrap(X,N)
    if (~is_octave())
        [phat,cint] = binofit(X,N);
    else
        [phat,cint] = binofit_dale(X,N);
    end
end
