% check if we are running octave
function b = is_octave();
    b = exist('OCTAVE_VERSION', 'builtin') ~= 0;
end
