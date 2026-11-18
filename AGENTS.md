# Work on these CFD examples

1. Read README, docs/BUILD.md, docs/MODELS_AND_LIMITS.md and the selected example README before changing code.
2. Keep upstream NGA2 at the recorded commit. Work in a separate case/build copy; preserve original inputs, failure logs and normal restart points.
3. Check existing solver ownership before launch. Never start a duplicate campaign. Never terminate unrelated processes. Use the local runner lock; coordinate external runs separately.
4. Start from a fresh input or verified normal checkpoint. Checkpoint compatibility requires the same source/binary identity, precision, partition and physical input; do not remove identity checks to load incompatible data.
5. Use real field values and native PLIC geometry for scientific results. Fix color ranges across a sequence. Keep explanatory illustrations separate from calculated results.
6. Do not clip velocity, temperature or resolved negative phase volumes to make a run pass. Find the stage/coordinate and repair consistent mass, momentum and energy coupling.
7. Retain pressure/projection, phase-volume, thermal and source conventions. Record divergence-minus-source, conservation, pressure residual, finite fields and phase positivity. MPI_ABORT may return zero.
8. Test one change at a time using short checks. Do not turn execution success into grid/time convergence or product accuracy. Known 2D velocity temporal accuracy is unresolved.
9. Keep source, run instructions and a concise transfer prompt. Do not commit raw runs, checkpoints, installed dependencies, private paths, account/SSH material or credentials.
10. Ask before pushing/publishing; local preparation is reversible. No obsolete HPC, memory or artificial file/read budgets are inherited from previous chats. Respect the new user's resource and usage limits.

Portfolio scope: also inspect examples/legacy_catalog.json and docs/LICENSE_MAP.md. The root LICENSE maps mixed terms; the NGA2 MIT notice does not relicense LICA Fortran. Complete build recipes and reference-only prototypes must remain distinct. Do not resume stopped hull work or old SSH controllers.

The documentedv3-to-v4 output-only restart is a narrow qualified exception: all8zero-step physical CSV files were byte-equal, followed by same-dt history/output checks. Do not generalize that exception to other binary/source changes.
