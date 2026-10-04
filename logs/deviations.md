# Deviations log

| Date | Item | Original | Changed to | Reason |
|---|---|---|---|---|
| 2026-10-03 | Boruta package version | latest CRAN (10.0.0) | 8.0.0, then **reverted to 10.0.0** the same day | v10 needs the Rust compiler (via `fru`). Pinned to 8.0.0 temporarily; cargo 1.93.1 was then installed and Boruta 10.0.0 installed. No deviation remains in effect. Recorded before any analysis ran. |
| 2026-10-03 | Adipose discovery cohort definition (GSE27951) | 11 ND / 10 T2D, n=21 (protocol/proposal) | GPL570 subseries, IGT excluded: 11 ND / 12 T2D, n=23 | GEO holds 33 human GPL570 samples (11 NGT / 10 IGT / 12 DM), not 21. IGT samples fit neither group. PI-level decision by the project owner (option 1). Group taken from the clinical-status field; 3 samples (GSM691135, GSM691139, GSM691140) have title vs field discordance, field used, awaiting confirmation from Dr. Thirumurugan. Proposal v3 and IMPLEMENTATION_PLAN counts (n=21) need updating. |
| 2026-10-05 | Methodology protocol | `docs/Methodology_Protocol_v1.md` referenced by CLAUDE.md but never written | Authored and locked 2026-10-05 (tag `protocol-v1`) | Written after GEO download and metadata inspection, before any expression data was normalised or analysed. Stated in the protocol header. |
| 2026-10-05 | Task 5 permutation self-test | 1 shuffle, accept AUC 0.40-0.60 | 20 shuffles; mean in 0.45-0.55, 95th pct < 0.65 | Single shuffle fails by chance ~19% (null SD ~0.075 at 20 vs 57). Changed before Task 5 was run. |
| 2026-10-05 | GSE27951 discordant labels | "awaiting PI confirmation" | Clinical-status field final; sensitivity rerun excluding GSM691135, GSM691140 | Decision by project lead. |
