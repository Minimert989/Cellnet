# Vendored CellNet core snapshot

The C++ implementation and public headers under `src/` are vendored so this
standalone R package builds without reaching into another checkout. The R
bridge is `src/bridge.cpp`; the core files are mirrored from the CellNet engine
repository and must be reviewed against that source when preparing a release.
The monorepo integration checkout carries the internal synchronization script;
it is intentionally not shipped as a working command in this standalone repo.
