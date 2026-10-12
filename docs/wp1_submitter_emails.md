# WP1: draft letters to the submitters and to GEO (NOT SENT)

Status: **DRAFT, not sent.** Written 2026-10-13 for the project lead (Juvas) to review, edit and sign. Supervisor review (Dr. Thirumurugan) is recommended before sending.
Evidence sources: `manifest/GSE38642_donor_overlap.csv` (from `scripts/02_donor_overlap.R`), reproduced independently by `scripts/16_wp1_overlap_map.R` (`results/stage8/16_overlap_pairwise.csv`, `results/stage8/16_overlap_donor_level.csv`).
Contacts come from the GEO series-matrix headers (`!Series_contact_*`), checked 2026-10-12:

| Series | GEO contact (header) | Institute (header) | E-mail in header |
|---|---|---|---|
| GSE41762 | Anders Rosengren | Lund University (CRC 91-11 SUS, Malmö, Sweden) | not listed; project lead to add a verified address |
| GSE38642 | Joao Fadista | Lund University Diabetes Centre (Malmö, Sweden) | joao.fadista@med.lu.se (as listed in the GEO header; check it is still current) |

Note before sending: the GSE38642/GSE41762 overlap has been noted in earlier literature (see project memory, novelty constraints). The letters therefore ask for confirmation and clarification only and do not present the overlap as a new discovery.

---

## Letter 1: to Dr Anders Rosengren (GSE41762)

**To:** Dr Anders Rosengren, Lund University (address to be completed)
**Cc:** Dr Joao Fadista (GSE38642 contact)
**Subject:** Question about the donors in GEO series GSE41762 and GSE38642

Dear Dr Rosengren,

We are a small computational group re-analysing public human islet expression data for type 2 diabetes. We use GSE41762, which you deposited, as our discovery cohort and are grateful that the data are openly available.

While checking the independence of public islet datasets, we compared GSE41762 with GSE38642, which comes from the same Nordic Islet Transplantation Programme material and from Lund University Diabetes Centre. We would be grateful if you could help us describe the two deposits correctly. What we observed:

1. **Shared arrays.** Each of the 63 CEL files in GSE38642 is byte-identical to one CEL file in GSE41762. We compared MD5 checksums of the decompressed CEL files. Each pair also has the same donor number in the sample title ("human islets, IDn"), so the 63 GSE38642 arrays seem to be a subset of the 77 GSE41762 arrays. The full list is in the appendix.
2. **Different diabetes labels for 11 donors.** For 11 of the 63 shared arrays, the diabetes status differs between the two series (table below). In 10 cases GSE41762 says "Diabetic donor" and GSE38642 says "non-diabetic donors"; in one case it is the other way round. For these donors, both deposits give the same age, BMI and HbA1c wherever HbA1c is present.
3. **Two smaller differences.** Donor ID30 (GSM1023682 / GSM946727) is recorded as female in GSE41762 and male in GSE38642. Donor ID48 (GSM1023700 / GSM946744) has HbA1c 5.6 % in GSE41762 and 5.9 % in GSE38642.

| Donor | GSE38642 sample | GSE41762 sample | Status in GSE38642 | Status in GSE41762 | HbA1c GSE38642 (%) | HbA1c GSE41762 (%) | Sex | Age | BMI |
|---|---|---|---|---|---|---|---|---|---|
| ID20 | GSM946717 | GSM1023672 | non-diabetic | diabetic | 6.1 | 6.1 | male | 57 | 23.9 |
| ID28 | GSM946725 | GSM1023680 | non-diabetic | diabetic | 6.2 | 6.2 | female | 59 | 36.6 |
| ID31 | GSM946728 | GSM1023683 | non-diabetic | diabetic | 6.2 | 6.2 | male | 59 | 26.3 |
| ID35 | GSM946732 | GSM1023687 | non-diabetic | diabetic | 7.0 | 7.0 | female | 64 | 28.4 |
| ID38 | GSM946735 | GSM1023690 | non-diabetic | diabetic | 8.6 | 8.6 | female | 54 | 27.0 |
| ID49 | GSM946745 | GSM1023701 | T2D | non-diabetic | not deposited | not deposited | male | 62 | 24.8 |
| ID54 | GSM946750 | GSM1023706 | non-diabetic | diabetic | 8.0 | 8.0 | male | 61 | 24.7 |
| ID56 | GSM946752 | GSM1023708 | non-diabetic | diabetic | 6.2 | 6.2 | female | 64 | 27.0 |
| ID61 | GSM946757 | GSM1023713 | non-diabetic | diabetic | 6.2 | 6.2 | male | 62 | 36.9 |
| ID62 | GSM946758 | GSM1023714 | non-diabetic | diabetic | 6.2 | 6.2 | female | 57 | 32.5 |
| ID63 | GSM946759 | GSM1023715 | non-diabetic | diabetic | 6.4 | 6.4 | male | 62 | 25.3 |

Our questions:

1. Can you confirm that the GSE38642 arrays are the same hybridisations as the corresponding GSE41762 arrays, from the same donors?
2. How was type 2 diabetes defined for GSE41762: by clinical diagnosis, by an HbA1c threshold (and if so, which one), or by both? The GEO record says "57 non-diabetic and 20 diabetic" but does not give the criterion.
3. For the 11 donors above, which status should a re-analysis use? Is there a reason, such as a different definition or updated clinical information, why they are labelled differently in the two series?
4. Would you support asking GEO to cross-reference the two series (for example with a note or a "related series" link), so that later users do not treat them as independent cohorts? We would only contact GEO with your agreement, or leave it to you if you prefer.

We are not suggesting that anything was done incorrectly. Different labels can follow from different, equally valid definitions. We want to describe the data accurately in our own work, and we will report your answer as you give it. If you prefer, we will simply say that we contacted you.

Our scripts are open, and we are happy to send them or a short report: `scripts/02_donor_overlap.R` and `scripts/16_wp1_overlap_map.R` (CEL checksums and covariate matching; no expression values are needed to reproduce the comparison).

Thank you for your time, and for making the data available.

Kind regards,

[Name, position — to be completed by the project lead]
[Institution, address]
[E-mail]
[Supervisor: Dr. Thirumurugan — confirm before sending]

---

## Letter 2: to Dr Joao Fadista (GSE38642)

**To:** joao.fadista@med.lu.se (from the GEO header; verify)
**Cc:** Dr Anders Rosengren (GSE41762 contact)
**Subject:** Question about the donors in GEO series GSE38642 and GSE41762

Dear Dr Fadista,

We are a small computational group re-analysing public human islet expression data for type 2 diabetes. Your series GSE38642 is among the most widely re-used islet datasets, and we are grateful that it is public.

While checking whether public islet datasets are independent of each other, we compared GSE38642 with GSE41762 (contact: Dr Anders Rosengren, Lund University). Both series come from the Nordic Islet Transplantation Programme. What we observed:

1. **Shared arrays.** Each of the 63 CEL files in GSE38642 is byte-identical (MD5 of the decompressed CEL file) to one of the 77 CEL files in GSE41762. Each pair also has the same donor number in the sample title ("human islets, IDn"). The full list is in the appendix.
2. **Different diabetes labels for 11 donors.** For 11 of the 63 shared arrays, the diabetes status differs between the two series. In 10 of them GSE38642 says "non-diabetic donors" and GSE41762 says "Diabetic donor"; donor ID49 is the other way round. Age, BMI and the deposited HbA1c are identical in both series for all 11 (table in Letter 1, repeated below).
3. **Two smaller differences.** Donor ID30 (GSM946727) is recorded as male in GSE38642 and female in GSE41762. Donor ID48 (GSM946744) has HbA1c 5.9 % in GSE38642 and 5.6 % in GSE41762.
4. **GSE50397 (also deposited by you).** 67 of the 89 CEL files in GSE50397 are byte-identical to GSE41762 arrays, with the same donor numbers in the titles. 56 of these are also GSE38642 arrays; the other 22 GSE50397 arrays (donors 75, 81, 82, 84–91 and 100–110) match no array in the other series we checked. In GSE50397 the deposited age is one year higher than in GSE41762 for 54 of the 67 shared donors, and two years higher for 2. BMI is the same as in GSE41762 for all 67, and HbA1c is the same wherever both series give it. Donor ID30 is recorded as male in GSE50397, as in GSE38642. GSE50397 has no diabetes-status field. We compared only file checksums and sample metadata for GSE50397 and did not open its expression data.

[Same 11-row table as in Letter 1]

Our questions:

1. Can you confirm that the GSE38642 arrays are the same hybridisations, from the same donors, as the corresponding GSE41762 arrays?
2. How was type 2 diabetes defined for GSE38642? The record says "54 non-diabetic and 9 diabetic". Several donors labelled non-diabetic have HbA1c of 6.5 % or higher (ID35 7.0 %, ID38 8.6 %, ID54 8.0 %). Was the label based on clinical diagnosis rather than HbA1c, or on a different threshold?
3. For the 11 donors above, which label should a re-analysis use?
4. For GSE50397: was age recorded at a different reference date (for example, calculated from year of birth), and which ND/T2D definition applies to its donors?
5. Would you support asking GEO to cross-reference GSE38642, GSE41762 and GSE50397, so that they are not treated as independent cohorts in later work? We would only contact GEO with your agreement, or leave it to you if you prefer.

We are not suggesting that anything was done incorrectly. We want to describe the data accurately and will report your answer as you give it.

We are happy to share our scripts (`scripts/02_donor_overlap.R`, `scripts/16_wp1_overlap_map.R`) and the per-sample comparison tables.

With thanks and kind regards,

[Name, position — to be completed by the project lead]
[Institution, address]
[E-mail]

---

## Note 3: to GEO curators (geo@ncbi.nlm.nih.gov) — send only after the submitters reply or after a reasonable wait

**Subject:** Possible cross-reference between GSE38642 and GSE41762 (and GSE50397) — shared CEL files

Dear GEO curators,

We would like to report an observation about two public series and ask whether a cross-reference would be appropriate.

- **GSE38642** (63 samples, GPL6244, contact Joao Fadista, Lund University Diabetes Centre) and **GSE41762** (77 samples, GPL6244, contact Anders Rosengren, Lund University) both describe human pancreatic islets from the Nordic Islet Transplantation Programme.
- All 63 CEL files in GSE38642 are byte-identical (MD5 of the decompressed CEL file) to 63 of the CEL files in GSE41762, one-to-one, and the paired samples carry the same donor number in their titles.
- For 11 of these pairs the deposited diabetes status differs between the two series, and for one pair the deposited sex differs. Sample accessions are listed below.
- [Optional, if the submitters agree] **GSE50397** (89 samples, GPL6244, same contact as GSE38642) contains 67 CEL files byte-identical to GSE41762 arrays (56 of them also in GSE38642). It is a SubSeries of GSE50398 and is not linked to GSE38642 or GSE41762 either.

The two series are not linked in GEO, so re-users may treat them as independent cohorts (for example, one as discovery and the other as validation). We have written to both submitters [date] to ask them to confirm and clarify. [Add their reply, or "no reply as of <date>".] We would be grateful if GEO could consider adding a note or a cross-reference between the two series. We leave any change to the sample annotations to the submitters.

Sample pairs with differing diabetes status (GSE38642 / GSE41762): GSM946717/GSM1023672, GSM946725/GSM1023680, GSM946728/GSM1023683, GSM946732/GSM1023687, GSM946735/GSM1023690, GSM946745/GSM1023701, GSM946750/GSM1023706, GSM946752/GSM1023708, GSM946757/GSM1023713, GSM946758/GSM1023714, GSM946759/GSM1023715.
Pair with differing sex: GSM946727/GSM1023682.
The full list of 63 pairs and our scripts are available on request.

Kind regards,

[Name, institution — to be completed by the project lead]

---

## Appendix: the 63 identical CEL pairs (GSE38642 → GSE41762)

Each pair carries the same donor ID in both series' sample titles (source: `manifest/GSE38642_donor_overlap.csv`, column `identical_CEL_GSE41762_gsm`).

| GSE38642 | GSE41762 | GSE38642 | GSE41762 | GSE38642 | GSE41762 |
|---|---|---|---|---|---|
| GSM946700 | GSM1023655 | GSM946701 | GSM1023656 | GSM946702 | GSM1023657 |
| GSM946703 | GSM1023658 | GSM946704 | GSM1023659 | GSM946705 | GSM1023660 |
| GSM946706 | GSM1023661 | GSM946707 | GSM1023662 | GSM946708 | GSM1023663 |
| GSM946709 | GSM1023664 | GSM946710 | GSM1023665 | GSM946711 | GSM1023666 |
| GSM946712 | GSM1023667 | GSM946713 | GSM1023668 | GSM946714 | GSM1023669 |
| GSM946715 | GSM1023670 | GSM946716 | GSM1023671 | GSM946717 | GSM1023672 |
| GSM946718 | GSM1023673 | GSM946719 | GSM1023674 | GSM946720 | GSM1023675 |
| GSM946721 | GSM1023676 | GSM946722 | GSM1023677 | GSM946723 | GSM1023678 |
| GSM946724 | GSM1023679 | GSM946725 | GSM1023680 | GSM946726 | GSM1023681 |
| GSM946727 | GSM1023682 | GSM946728 | GSM1023683 | GSM946729 | GSM1023684 |
| GSM946730 | GSM1023685 | GSM946731 | GSM1023686 | GSM946732 | GSM1023687 |
| GSM946733 | GSM1023688 | GSM946734 | GSM1023689 | GSM946735 | GSM1023690 |
| GSM946736 | GSM1023691 | GSM946737 | GSM1023692 | GSM946738 | GSM1023693 |
| GSM946739 | GSM1023694 | GSM946740 | GSM1023695 | GSM946741 | GSM1023696 |
| GSM946742 | GSM1023698 | GSM946743 | GSM1023699 | GSM946744 | GSM1023700 |
| GSM946745 | GSM1023701 | GSM946746 | GSM1023702 | GSM946747 | GSM1023703 |
| GSM946748 | GSM1023704 | GSM946749 | GSM1023705 | GSM946750 | GSM1023706 |
| GSM946751 | GSM1023707 | GSM946752 | GSM1023708 | GSM946753 | GSM1023709 |
| GSM946754 | GSM1023710 | GSM946755 | GSM1023711 | GSM946756 | GSM1023712 |
| GSM946757 | GSM1023713 | GSM946758 | GSM1023714 | GSM946759 | GSM1023715 |
| GSM946760 | GSM1023716 | GSM946761 | GSM1023717 | GSM946762 | GSM1023718 |
