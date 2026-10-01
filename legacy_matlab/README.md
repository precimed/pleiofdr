# Legacy MATLAB implementation

The original MATLAB/Octave implementation of pleioFDR, kept for reference. It is no longer
maintained; use the Python package (`pleiofdr --config config.txt`, see the main README).

It still runs from the repository root, where the config files and data live:

```
matlab -nodisplay -nosplash -r "addpath('legacy_matlab'); runme; exit"
```

`tests/fixtures/make_golden.m` uses it to produce the reference outputs that the Python test suite
compares against. MIGRATION_NOTES.md maps each MATLAB function to its Python counterpart.
