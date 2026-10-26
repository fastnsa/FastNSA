# Troubleshooting

## Native Library Not Found

Symptom: import reports that libfastnsa.so, libfastnsa.dylib, or fastnsa.dll
cannot be found.

1. Confirm that pip installed a wheel compatible with the active interpreter.
2. Print the interpreter and package locations.
3. Reinstall in the same environment.
4. For a developer build only, set FASTNSA_LIBRARY_PATH to the exact library.

~~~bash
python -c "import sys; print(sys.executable); print(sys.version)"
python -m pip show fastnsa
python -m fastnsa._selftest
~~~

## Library Exists but Cannot Be Loaded

An OSError usually indicates a platform, architecture, or dynamic dependency
mismatch. Check that Python and the library are both 64-bit and were built for
the same OS and CPU architecture. On Windows, use the correct MSVC runtime. On
Linux and macOS, inspect the library dependencies with platform tools.

Do not rename a library from another platform to the expected filename.

## Optional Package Import Error

Install the extra named in the error message using the same Python executable:

~~~bash
python -m pip install "fastnsa[geo]"
~~~

Inside QGIS, run installation against QGIS's Python environment. A successful
installation in a Conda or system shell may not be visible to QGIS.

## Distances Are Implausible

The most common cause is a geographic CRS. Coordinates in EPSG:4326 are degrees,
not metres. Project both network and events to an appropriate projected CRS,
rebuild the network, and repeat mapping. Check that edge lengths and all
distance parameters use the same unit.

## Mapped Edge IDs No Longer Match

Mapped edge IDs depend on native edge order. They become invalid after edge
sorting, filtering, topology edits, reprojection/reconstruction, or switching
to a different network. Store and compare a stable network signature before
reusing mappings.

## ST K-function Rejects Times

Check the following:

- Every event has a finite numeric time.
- At least two events exist.
- The maximum and minimum times differ.
- Temporal thresholds use the same unit as event times.
- MTS and EARTH receive at least two evenly spaced thresholds per axis.

Datetime columns loaded through PointEvents are converted to Unix seconds. If
you want thresholds in days, explicitly convert the event array to days or use
second-valued thresholds consistently.

## No Theoretical Curve or Envelope

This is expected when simulations is zero. FastNSA estimates k_theo from Monte
Carlo simulations; it does not insert an analytic curve. Set simulations to a
positive count.

## Cancellation Is Delayed

Cancellation is cooperative. The backend checks callbacks at safe boundaries,
so a long native stage may finish its current unit of work before returning.
Do not forcibly unload the native library or destroy the network while a task
is running.

## Memory Use Is High

- Set return_simulations=False unless every simulation is required.
- Reduce ST spatial or temporal grid resolution.
- Increase lixel_length to create fewer NKDE segments.
- Process exploratory runs with fewer simulations, then increase for final
  inference.
- Reuse one network and one mapped point set instead of rebuilding them.

## Analysis Is Slow

Measure stages separately before optimizing. Typical costs include graph or
coordinate conversion, nearest-edge mapping, lixel generation, threshold evaluation, and
Monte Carlo repetition. Use projected, valid, reasonably simplified geometry;
reuse precomputed mappings; compare optimized methods on representative data;
and avoid threshold ranges far beyond the study scale.

Release wheels should use optimized native builds. Local Debug builds can be
substantially slower.

## Native Crash or Signal 11

A Python exception should normally be raised for invalid public API inputs. A
segmentation fault is a defect or binary mismatch, not a normal validation
outcome.

Record the FastNSA version, OS, CPU, Python version, method, parameters, point
count, edge count, whether simulations or weights were enabled, and a minimal
reproducible dataset. Re-run with the current release and self-test before
filing an issue. Avoid continuing to use the same process after native memory
corruption is suspected.

