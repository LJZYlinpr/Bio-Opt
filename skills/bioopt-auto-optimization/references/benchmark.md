# Paired benchmarks and evidence

Run in an allocated Linux compute environment. Prepare two noninteractive commands that write separate outputs and an independent validation command. Relative paths resolve from the current directory. A thin adapter can wrap an official entry point. Archive the adapter and its arguments.

```json
{
  "case_id": "algorithm-benchmark-example",
  "claim_scope": "declared workflow including setup and output writing",
  "equivalence_contract": "all declared output fields match exactly",
  "inputs": ["frozen-input.npz"],
  "baseline": ["python", "run_reference.py", "frozen-input.npz", "{output}/result.npz"],
  "candidate": ["python", "run_candidate.py", "frozen-input.npz", "{output}/result.npz"],
  "validate": ["python", "/path/to/skill/scripts/compare_numpy.py", "{baseline}/result.npz", "{candidate}/result.npz"],
  "threads": 1,
  "timeout_seconds": 300,
  "warmups": 1,
  "repeats": 3
}
```

```bash
python /path/to/skill/scripts/paired_commands.py case.json --out evidence/pair-01
```

To measure baseline repeatability, give the [paired command runner](../scripts/paired_commands.py) two identical baseline commands. Use a fresh directory for formal comparisons and validate each pair of outputs independently. The runner hashes inputs, which can warm the file cache. For cold-start comparisons or inputs shared across runs, specify the warm-up policy and report any measured shared setup costs and the resulting workflow time separately. Peak RSS comes from Linux `wait4`. Measure device peaks and aggregate memory across processes with workload-specific tools. Command timing includes startup and writing outputs. Validation is timed separately.

The [NumPy comparator](../scripts/compare_numpy.py) defaults to comparing output keys, shapes, dtypes, and bit patterns. It rejects nonfinite floating-point values and pickled data. Use `--float-key` with `--atol` or `--rtol` only for fields with explicitly approved tolerances. For output formats with special semantics, specify the acceptance criteria and use an independent format-specific comparator. A validation failure must return a nonzero exit status.

## Records to save

- The upstream example, commit, input hashes, commands, environment, patch, device, and resources.
- The timing boundary, included setup, run order, raw measurements, median and variability, and peak memory.
- The frozen scientific criteria, independently checked fields and differences, rejected candidates, and grounds for acceptance.
- For batch runs, save task, implementation, and condition records. Use explicit completion records to resume runs and preserve original errors. If resources permit one pair only, report the count and treat the result as provisional.
- For a comparison with ordinary Codex, use the same starting commit, scientific criteria, compute, time, and candidate trial budget on both sides. Apply the same frozen verifier, then compare the best accepted results.
