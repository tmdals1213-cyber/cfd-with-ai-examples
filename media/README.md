# Shareable media

- [3D landscape, 24 s, 1920×1080](IMMERSION_COOLING_LANDSCAPE_DRAFT01.mp4)
- [3D portrait, 24 s, 1080×1920](IMMERSION_COOLING_PORTRAIT_DRAFT01.mp4)
- [Example overview, 20 s](EXAMPLE_OVERVIEW.mp4)

The 3D draft uses twelve actual native XYZ PLIC snapshots from 1–20 ms and the actual final-temperature/velocity slice. The English synthetic voice is generated at a brisk native speech rate, with no time-stretch. The full MP4s were decoded for verification. The overview is a slide summary of recorded results, not a new time-resolved CFD animation.

[Storyboard](STORYBOARD.md) · [Narration](NARRATION_EN.txt) · [Data provenance](../provenance/video_draft01.json) · [Scope](../docs/MODELS_AND_LIMITS.md)

Re-render requires the original complete run files (intentionally not bundled). Example command:

```bash
python3 -m pip install -r requirements-render.txt
python3 scripts/render_showcase.py --run /path/to/completed_case --execution /path/to/EXECUTION.json --limit .02 --output /path/to/new_video --ffmpeg ffmpeg
# Optional local Windows narration, before rendering with --audio-dir:
# powershell -File scripts/voice_windows.ps1 -OutputDirectory /path/to/new_audio
```

Changing --limit alone does not synthesize a longer result; complete matching native PLIC, history and all eight field packets are required. A later 100 ms result must receive a new manifest, not a relabelled old movie. The supplied render layout assumes this 80³ geometry and case; adapt it deliberately for other examples.

For the existing preview gallery alone: `python3 scripts/render_overview.py --ffmpeg ffmpeg --output /path/to/new_overview.mp4`. This does not require raw fields.

## Earlier project families

[Ten-car landscape](VEHICLE_CAR10_LANDSCAPE.mp4) · [ten-car portrait](VEHICLE_CAR10_PORTRAIT.mp4) · [EV](VEHICLE_EV_PAIR.mp4) · [SUV](VEHICLE_SUV_PAIR.mp4) · [DrivAer wake](VEHICLE_DRIVAER_WAKE.mp4). The ten-car English sequence is a camera orbit around one matched final field, not a new time history.

[Runner landscape](RUNNER_RECORDED_SUMMARY_LANDSCAPE.mp4) · [runner portrait](RUNNER_RECORDED_SUMMARY_PORTRAIT.mp4). These show recorded offsets and historical pressure-component means. They do not recover the deleted old flow footage or establish total drag reduction.
