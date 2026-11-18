# Continuous saved-field playback

`render_showcase.py` is the existing narrated20ms introduction. `render_continuous_flow.py` prepares a continuous result view: native3D PLIC, a temperature slice and velocity magnitude/directions at the same saved time. It preserves the45–63C and0–0.6m/s display ranges. It repeats actual saved states in slow motion and does not interpolate CFD fields or invent intermediate steps.

```bash
python3 -B scripts/render_continuous_flow.py --run /path/to/accepted-run --output /path/to/new-video --ffmpeg /path/to/ffmpeg --limit 0.1 --seconds 12 --form landscape
```

Use `--form portrait` for the vertical layout. Install requirements-render.txt in a separate environment. The current reader deliberately checks the80cubed/20mm geometry; adapt it explicitly for a different case.

Before rendering, all eight physical CSV packets at one matched time must byte-match exported float32temperature, latent source, VOF and velocity. Accepted history residuals, phase ranges and the independent gas-volume sum are checked. Native PLIC and field schedules must agree. The local25–35ms five-state QA passed these checks and a complete video decode. This is a renderer/data-integrity check, not a claim that100ms or a complete cooling cycle has finished.

The continuation target is first100ms with2.5ms saved intervals. Inspect actual phase geometry and motion before deciding whether a longer300ms stage is useful. If the solver stops, preserve the failure and resume only a qualified normal state after diagnosis.
