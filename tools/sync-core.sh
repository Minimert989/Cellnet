#!/bin/sh
set -eu

PACKAGE_DIR=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
REPOSITORY_ROOT=$(CDPATH= cd -- "$PACKAGE_DIR/../.." && pwd)

for file in binary.cpp compiler.cpp engine.cpp parser.cpp scc.cpp synthetic.cpp temporal.cpp; do
    cp "$REPOSITORY_ROOT/src/$file" "$PACKAGE_DIR/src/$file"
done

mkdir -p "$PACKAGE_DIR/src/cellnet"
for header in "$REPOSITORY_ROOT"/include/cellnet/*.hpp; do
    cp "$header" "$PACKAGE_DIR/src/cellnet/$(basename "$header")"
done

printf '%s\n' "Synchronized CellNet C++ source and headers from $REPOSITORY_ROOT"
