# K-band M/L profile fits

This utility reads ML_Andreea.csv and fits each (Sample, Galaxy) group with constant, linear, and log-linear unweighted profiles. Linear models are centered on the group mean Radius_kpc; invalid or non-positive M/L rows are recorded as exclusions.

Run the production calculation:

    python3 fit_ml_profiles.py --input ML_Andreea.csv --output-dir outputs/ml_profile_fits

Run the dependency-free validation:

    python3 test_ml_profiles.py

The test uses only Python's standard library plus the runtime pandas dependency already required by the fitting script. It checks known constant, linear, and log-linear profiles and verifies invalid-row accounting. Outputs are ml_profile_fits.csv, ml_profile_exclusions.csv, ml_profile_summary.csv, and ml_profile_metadata.json. COSMOS_kin remains identifiable in per-group records while COSMOS summaries combine COSMOS_DR1 and COSMOS_kin. The utility does not calculate V_disk or RAR quantities.
