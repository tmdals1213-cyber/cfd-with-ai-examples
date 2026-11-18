# Vehicle / runner source transfer

## Vehicle and cylinder Fortran

The source keeps the original LICA copyright/noncommercial notices. Install GNU Fortran with OpenMP. Use `python3 scripts/legacy_case.py list`, then `build` one selected complete Fortran case. The build uses legacy default8-byte integer/real conventions and the explicitly ordered modules. It is isolated under ignored build/. Windows adds the original stack-linker option; Linux uses its normal linker.

`prepare` copies only retained small input files into a new ignored runs/ directory and changes IREAD=0/IRESET=1 for a new cold start. It does not run a solver, import a failed state, modify the originals or fetch any geometry. Preserve the physical/nondimensional units and review the remaining input values before running.

Geometry and IBM preprocessing must be recreated for the exact body/grid. The retained preprocessor reference is specialized at compile time (161×81×81, LES/no interpolation); adapting dimensions, LES and interpolation flags is required for another variant. It expects an appropriately generated body.sdf; it does not magically convert a source STL into a certified watertight union. The original artist models, mesh repair and orientation/scale records are external and attributed in VEHICLE_ASSET_LEDGER.json. A preprocessor binary or SDF from another grid is not a reusable operand.

Legacy source may contain platform/runtime assumptions beyond compilation. A successful build is not a fresh simulation, a tested cold start, or physical validation. Before running a legacy case, coordinate it with the active immersion calculation and inspect the small setup/asset requirements. No old SSH/account route is supplied.

## Final runner C++

Build `15_runner_final` with GNU C++ / OpenMP. A non-Windows symbol-export macro was added; the mathematical kernel is retained. The cold5×7×9 and preview76×29×72 libraries are built separately. Small ABI/analytic Vreman checks advance zero fluid steps.

The original reference_pipeline is an exact-artifact stage interface with configuration/geometry/motion dependencies. It is included for agent integration, not as a new turnkey launch command. Its old per-stage control/once policies are historical implementation details, not inherited authorizations. Adapt it around a fresh configuration and the external motion/topology/geometry provider; do not execute former remote controllers. requirements-legacy.txt lists the Python mathematical libraries.

The normal43-sample motion NPZ and topology/geometry source are available only in the user's local private original archive. No installed DLL, full fluid array or raw restart is embedded. The private local migration map is separate from the public candidate.

## Earlier / stopped work

The four-runner prototype, early chip thermal, acoustic/hairdryer preparation and stopped hull/VOF sources are explicit reference-only entries. Their inclusion does not assert successful execution or restart authorization. The hull work stays stopped. All current phase-change cooling work remains under the separate NGA2 case.
