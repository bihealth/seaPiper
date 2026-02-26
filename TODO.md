# TODO

## Next Tests (Priority Order)

1. `Rseasnap` constructor integration  
Target `R/make_seapiperdata.R` (`seapiperdata_from_rseasnap`) and `R/data_prep.R` (`.prepare_data_single_pipeline`).  
Add mocked `get_*` tests for happy path and partial failures (`tmod_*` unavailable, missing `dataset_title`, contrast column filtering).

2. `rld` extraction edge cases  
Target `R/data_prep.R` (`.extract_rld_matrix`).  
Add S4 assay-shape fixtures and failure cases to verify exact error behavior.

3. Multi-dataset YAML and path resolution  
Target `R/make_seapiperdata.R` (`seapiperdata_from_yaml`) and path helpers.  
Test multiple datasets, absolute vs relative paths, and mixed literal/object YAML values.

4. Merge/schema edge cases  
Target `R/make_seapiperdata.R` (`merge_seapiperdata`) and normalization helpers.  
Test asymmetric fields, all-`NULL` field collapse, and malformed named-list fields.

5. Remaining server wrappers with reactivity  
Target `R/server_gene.R`, `R/server_tmod.R`, `R/server_info.R`.  
Use `shiny::testServer()` to verify `observeEvent` behavior and output rendering on input changes.

6. UI fallback behavior  
Target `R/ui.R` (`.pipeline_dashboard_body`).  
Test fallback generation of `pipeline_choices` and `cntr_titles` when config/titles are missing.

7. End-to-end app smoke tests  
Add `shinytest2` smoke tests: app boots with minimal/full data, expected tabs appear/disappear, and no runtime errors.

8. Deprecated alias coverage  
Target deprecated aliases in `R/make_seapiperdata.R` (`make_seapiperdata`, `rseasnap_to_seapiperdata`, `custom_yaml_to_seapiperdata`).  
Assert aliases call new functions and emit deprecation warnings.
