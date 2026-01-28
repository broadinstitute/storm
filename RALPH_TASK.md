---
task: Build STORM (Structural & Tandem-Repeat Optimized Regression Models)
test_command: "cargo test"
---

# Task: STORM

Build STORM, a Rust crate with a Python front-end for association testing of
structural variants and tandem repeats using long-read data.

---

## Success Criteria

1. [x] Integrated SV VCF can be parsed (SVTYPE, SVLEN, GT)
2. [x] TRGT VCF can be parsed (TRID, AL, GT)
3. [x] TRExplorer BED file can be ingested
4. [x] TRExplorer JSON file can be ingested
5. [x] TRExplorer BED and JSON can be joined into a unified catalog
6. [x] SV records can be mapped to catalog repeat loci by overlap
7. [x] TestUnit objects can be created for SVs
8. [x] TestUnit objects can be created for repeat-proxy loci
9. [x] TestUnit objects can be created for true repeats
10. [x] Resolver separates presence from allele values
11. [x] Repeat-proxy SV alleles can be grouped by locus
12. [x] Diploid repeat lengths can be reconstructed from INS and DEL alleles
13. [x] TRGT allele lengths override proxy alleles when available
14. [x] Resolved genotypes record presence source
15. [x] Resolved genotypes record allele source
16. [x] Arrow cache can be written
17. [x] Parquet cache files can be written
18. [x] Cache includes test units table
19. [x] Cache includes genotypes table
20. [x] Cache includes catalog table
21. [x] Cache includes features table
22. [x] Cache includes provenance metadata
23. [x] storm explain prints resolved genotype details
24. [x] Plan YAML can be parsed
25. [x] Plan rules deterministically select encodings
26. [x] Plan rules deterministically select models
27. [x] S encoding (L1+L2) is implemented
28. [x] M encoding (max allele) is implemented
29. [x] D encoding (|L1-L2|) is implemented
30. [x] Tail indicator encoding is implemented
31. [x] Categorical bin encoding is implemented
32. [x] Internal StormGLM backend is implemented
33. [x] Linear regression is supported
34. [x] Logistic regression is supported
35. [x] Categorical regression is supported
36. [x] BinomiRare test is supported
37. [x] Firth logistic regression is supported or feature-flagged
38. [ ] Covariates including PCs are supported
39. [ ] Rare-variant ladder logic is implemented
40. [ ] Association results can be written to Parquet
41. [ ] Results include statistics and p-values
42. [ ] Results include counts and call rates
43. [ ] Results include encoding and model metadata
44. [ ] Python API can build cache
45. [ ] Python API can load cache via Polars
46. [ ] Python API can run StormGLM
47. [ ] Python API can call explain
48. [ ] Jupyter notebook demonstrates end-to-end run
49. [ ] All tests pass
50. [ ] .ralph/progress.md contains DONE

---

## Context

- Primary input is an integrated long-read SV VCF
- Overlay input is a TRGT VCF for repeat allele lengths
- Repeat catalog is TRExplorer BED plus JSON
- Storage and interchange use Arrow and Parquet
- Rust core prioritizes correctness
- Python layer uses Polars for interactivity
- Mixed models are explicitly out of scope
- Design must allow future HLA and KIR overlays