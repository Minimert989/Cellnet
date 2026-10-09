# Vendored CellNet core snapshot

The C++ implementation and public headers under `src/` are a source snapshot
from this monorepo. They are included so `R CMD build` can produce a
self-contained R package tarball. Refresh them with:

```sh
sh integrations/cellnetR/tools/sync-core.sh
```

Review the copied C++ diff before release. The R bridge is `src/bridge.cpp`;
the vendored core files should otherwise match the repository's `src/` and
`include/cellnet/` files. No solver changes belong in this integration package.
