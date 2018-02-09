# Number of samples
total_samples = 1477
eff_samples = 1449 #no kidney
eff_samples = 1339 #no kidney liver

# Prefiltering: delete splicing events (SE) with >=0.1 NA's
* **psi.tissue**: Merge PSI + description (by sampleId) *NOT in GITHUB

* **notNA_tissue.tsv**: percentage of NA per SE and tissue (not kidney) *with dcast()
-> delete NA if percentage >= 0.1

* **depurated_notNA_tissue.tsv**: percentage of NA per SE and tissue, after deleting SE rows that percentage >=0.1 (only keeping SE <=0.1 NAs).

* **psi.tissue.filtered** (sample_id, splicing_event, PSI, tissue): filtered of NAs and no kidney *NOT in GITHUB

(*) graphics: 
- notNA_tissue vs. depurated_notNA_tissue
- psi.tissue vs. psi.tissue.filtered

# Training and Testing sets
training_set = 100/tissue * 5 = 34.50656% 

* **training_set** = 100 samples/tissue -> 500 samples *NOT in GITHUB
> (500/1449)*100
34.50656

* **testing_set** = 1449-500 = 949 *NOT in GITHUB
> (900/1449)*100
62.1118

