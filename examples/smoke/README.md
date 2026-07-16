# SCP-ResamplerIO smoke example

This no-argument example installs the Release tool and validates the published
Hemisphere Re1370 sample through `scp_nek5000_smoke`. The sample is the middle
20 Nek5000 source frames (490 through 509), resampled over the full source
bounds to a `250 x 80 x 50` structured grid.

The validator reads all 20 VTI frames and checks geometry,
`Velocity(float3)`, the binary valid-point mask, time ordering, PVD references,
and provenance JSON. The prepared data has 993,170 valid samples out of
1,000,000 in every frame.

Windows:

```powershell
.\examples\smoke\run.ps1
```

Linux:

```bash
bash ./examples/smoke/run.sh
```

The private `hweifluids/SCP-examples` repository versions the dataset under
`hemi_resampled/`; only its 20 VTI payloads use Git LFS. In the Streamcenter+
superproject the runner uses `examples/SCP-examples`. A detached
`SCP-ResamplerIO` checkout uses a sibling `SCP-examples` checkout when present.
If neither exists, authenticated GitHub CLI downloads and verifies the private
`smoke-data-v1` Release ZIP.

Prepared test data is copied to the fixed ignored path
`.tests/smoke/SCP-ResamplerIO/data` at the superproject root, or to the same
relative path below the detached module. No dependency, data, build, or output
path parameter is accepted, and no project environment variable is required.

The original 1,000-frame Nek5000 sequence is maintainer-only source provenance;
it is not needed by this portable smoke test. The `scp_nek5000_smoke` generator
remains available for reproducibly preparing a future version of the example
asset, while the generic `Resample_solver` continues to accept normal
`.nek5000` operation decks without being limited to these frame indices or
dimensions.
