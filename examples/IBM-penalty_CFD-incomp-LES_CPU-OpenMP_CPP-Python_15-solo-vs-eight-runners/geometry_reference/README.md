# Moving union geometry provider reference

Retained source for the original runner geometry provider, plus nearest-triangle and expansion-predicate headers. The binary and third-party topology/motion inputs are excluded. This source provides winding/union queries, segment events and moving-wall geometry; it is not a flow solver.

Windows C++17 shared-library build (MinGW):

```bash
g++ -std=c++17 -O2 -shared rnu_provider_local_utf8.cpp -o original_union.dll
```

For a Linux port, replace the Windows export annotation with default-visibility C ABI and verify int32/float64/uint64/EdgeEvent layouts. Inspect the `rnu1_rd_*` exports; supply the retained topology, pose/motion ABI and exact epoch/stage times explicitly. Do not substitute another runner mesh without rebuilding/adapting the inputs.

Included as source provenance/reference. The portfolio build matrix does not claim this provider or a full moving-geometry pipeline was rebuilt. No geometry/flow body was executed during packaging.
