# Reproduction limits and required private inputs

- The agreed JSD implementation is the supplied range-6 transfer package. One real MRI was tested end to end; no full-cohort FreeSurfer-to-JSD recomputation was performed.
- The complete 25-checkpoint nested classifier training was not repeated. Fold generation, saved training/validation preprocessing, held-out probabilities and source function equivalence were checked. OOF evaluation and fixed scoring were executed.
- The complete connection fitting batch was not repeated. One connection was refitted with both formal models; all saved formal connection fits were collected and their full inference families compared.
- Full execution requires the authorized reconstructed scan directories, a reviewed local manifest, the curated MRI-centered clinical table and sufficient local computation. These inputs are not missing code and are not public repository content.
- The MRI-centered clinical table is an input boundary. No raw clinical database reconstruction, new assessment matching rule or new date window is supplied.
- The manuscript's surface display reference is restricted and is not bundled. To reproduce its exact panel A appearance, use that authorized local reference. Another DKT surface gives a different rendering of the same regional statistics.
- Conda environment creation, non-Windows platforms and different numerical-library versions were not tested. The actual installed dependencies loaded successfully; exact tested versions are recorded.
- The research team must assign authorship and a software license before public code distribution. Packaging does not grant rights to redistribute ADNI data or fitted model artifacts.
