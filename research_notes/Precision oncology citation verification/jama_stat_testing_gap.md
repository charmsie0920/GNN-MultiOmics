# Verification of JAMA Network Open NGS-disparities paper and its STAT News coverage

## Does DOI 10.1001/jamanetworkopen.2026.5585 resolve to the claimed paper, and is it real?

### Takeaway
The DOI is real and resolves correctly. It is NOT broken. It points to: Hage Chehade C, Ozay ZI, Jo Y, et al. "Trends and Disparities in the Use of Next-Generation Sequencing in Patients With Cancer in the United States." *JAMA Netw Open.* 2026;9(4):e265585. Published online April 7, 2026 — a Flatiron Health EHR cohort of 63,294 patients with advanced/metastatic breast, prostate, NSCLC, colorectal, and pancreatic cancers, exactly matching the claimed description. A correction was issued May 7, 2026.

### Cited Findings
- `https://doi.org/10.1001/jamanetworkopen.2026.5585` returns an HTTP 302 redirect to `https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397`, which is a live JAMA Network Open article page — not a 404 — [DOI resolver / JAMA Network Open](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- Title on the live page: "Trends and Disparities in the Use of Next-Generation Sequencing in Patients With Cancer in the United States" — [JAMA Network Open](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- Authors: Chadi Hage Chehade, Zeynep Irem Ozay, Yeonjung Jo, et al. (17 authors total), senior/last author Umang Swami — [JAMA Network Open](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- Full citation: *JAMA Netw Open.* 2026;9(4):e265585, published online April 7, 2026 — [JAMA Network Open](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- Cohort size as stated on the page: "63 294 patients with advanced or metastatic cancer were eligible and included," drawn from 280 US cancer clinics via the Flatiron Health EHR-derived database, diagnosed 2018–2022 — [JAMA Network Open](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- A correction notice (dated May 7, 2026) amended wording in the Abstract and Conclusions, changing "most patients did not undergo tumor genomic testing" to "many patients did not undergo" tumor genomic testing — [search-result summary of JAMA correction notice](https://www.managedhealthcareexecutive.com/view/most-patients-with-advanced-cancer-still-not-receiving-genomic-testing-study-finds); independently corroborated by a second summary citing the same correction language — [search aggregation of JAMA correction text]
- Independent secondary confirmation that this is a real, widely covered study: Targeted Oncology (April 14, 2026), AJMC (April 7, 2026), UroToday abstract page, EMJ Reviews, Medscape, Healio (April 22, 2026), Managed Healthcare Executive (April 17, 2026), KFF Health News, and Georgetown press materials all describe the same Chehade et al. JAMA Network Open study with matching cohort size and topic — [Targeted Oncology](https://www.targetedonc.com/view/ngs-uptake-rises-but-equity-lags-behind-across-us-cancer-clinics); [AJMC](https://www.ajmc.com/view/sdoh-drive-disparities-in-ngs-access-across-cancers); [Medscape](https://www.medscape.com/viewarticle/tumor-genomic-testing-too-little-too-late-2026a1000b7g); [Managed Healthcare Executive](https://www.managedhealthcareexecutive.com/view/most-patients-with-advanced-cancer-still-not-receiving-genomic-testing-study-finds)

### Inferences
- The DOI, title, journal, year, and cohort size in the claimed citation are all accurate as given. There is no indication this is a fabricated or mismatched DOI.
- The May 2026 correction (softening "most" to "many") is itself evidence that the overall non-testing share is closer to "about half" or "a plurality" rather than comfortably over 50%, which is directly relevant to evaluating the "over 40%" claim below.

### Gaps
- I was not able to independently retrieve a PMID/PMCID for this article (PubMed- and PMC-site searches for the DOI and lead author did not surface a direct record in the time available). The JAMA Network Open publisher page itself is nonetheless a solid primary-source confirmation, so this gap does not undermine the authenticity finding, but a presenter wanting a PMC link should look it up directly on pubmed.ncbi.nlm.nih.gov.

## Does the paper report that "over 40%" (or roughly half) of patients received no genomic/NGS testing at all? What is the precise figure and wording?

### Takeaway
The paper does **not** state a single headline "over 40%" (or "half") figure for the whole cohort in those terms. Its reported numbers are broken out **by cancer type**, and they range from 33.0% (colorectal) to 59.7% (metastatic breast) with no testing at all. A simple weighted average across the five cancer types (computed from the per-type counts, not stated by the authors themselves) comes to approximately 41.3% of the full 63,294-patient cohort receiving no testing — which does support "over 40%" as an aggregate, and is independently consistent with secondary press coverage, but it is a derived number, not a sentence printed in the paper. The paper's own qualitative language (post-correction) is "many patients ... did not undergo tumor genomic testing."

### Cited Findings
- Per-cancer-type "no testing" figures as reported in the article's Table 1 (via direct retrieval of the JAMA Network Open full-text page): metastatic breast 7,215/12,119 = 59.7% untested; metastatic prostate 2,179/4,341 = 50.2%; advanced NSCLC 9,221/27,034 = 34.1%; metastatic colorectal 4,506/13,654 = 33.0%; metastatic pancreatic 3,008/6,166 = 48.8% (denominators back-calculated from percentages/counts as extracted) — [JAMA Network Open](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- Summing the "no testing" counts across all five cancer types and dividing by the full 63,294-patient cohort gives ≈41.3% (26,129/63,294) — this is an arithmetic combination performed during retrieval, not a sentence quoted verbatim from the paper — [derived from JAMA Network Open Table 1 data](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- Independent secondary corroboration of essentially the same aggregate figure: Medscape's coverage states "found that 41.3% had no genomic testing" for the ~63,000-patient cohort — [Medscape](https://www.medscape.com/viewarticle/tumor-genomic-testing-too-little-too-late-2026a1000b7g)
- The paper's own post-correction qualitative wording is "many patients with advanced or metastatic cancers did not undergo tumor genomic testing" (changed from the original "most patients," per the May 7, 2026 correction) — [JAMA Network Open correction, as summarized in secondary coverage](https://www.managedhealthcareexecutive.com/view/most-patients-with-advanced-cancer-still-not-receiving-genomic-testing-study-finds)
- Other secondary outlets describe the gap in varying, sometimes inconsistent framings: Healio (April 22, 2026) says "nearly half" of patients did not undergo NGS; Inside Precision Medicine describes "40% to 50% of patients still weren't getting testing" by the final study year — [Healio]; [Inside Precision Medicine](https://www.insideprecisionmedicine.com/topics/precision-medicine/u-s-tumor-testing-suffers-under-access-disparities/)

### Inferences
- "Over 40%" is a defensible, roughly accurate characterization of the pooled non-testing rate across the five cancer types in this cohort, and it is corroborated by at least one independent secondary source doing the same arithmetic (Medscape's 41.3%). However, it should not be presented as a verbatim quotation from the JAMA paper — the paper itself does not appear to print one single aggregate percentage in that form; it reports per-cancer-type rates, and its own prose says "many patients," not a specific percentage.
- For a formal academic citation, the more defensible and precise statement is either (a) cite the per-cancer-type range (33.0%–59.7% untested) directly from the paper, or (b) cite the ~41.3% pooled figure as a derived/calculated statistic (optionally citing Medscape as having independently calculated the same number), rather than attributing a single "over 40%" sentence to the paper's own text.

### Gaps
- I did not locate the exact wording of the paper's Abstract/Key Points sentence verbatim (only a retrieval-tool paraphrase/extraction of the full-text page); a presenter should pull the exact Abstract sentence directly from the JAMA Network Open page (and note the May 7, 2026 correction's wording change) before quoting it verbatim in a slide.

## Does the paper report median time-to-testing by cancer type, with NSCLC fastest and metastatic prostate cancer slowest at ~11.0 months? What is the EXACT NSCLC figure — is "41 days" supported?

### Takeaway
**"41 days" is NOT supported by this paper and should not be used.** The article (via direct retrieval of the full-text page) reports the median time to NGS for advanced NSCLC as **0.96 months (95% CI, 0.93–0.96 months)** — which converts to approximately **29 days**, not 41. The claim that NSCLC is fastest and metastatic prostate cancer is slowest at roughly 11.0 months is correct: metastatic prostate cancer's median time to NGS is reported as **11.0 months (95% CI, 9.8–11.0 months)**, the longest of the five cancer types, and NSCLC's 0.96 months is the shortest. So "under 1 month" is an accurate qualitative description for NSCLC, but the specific figure "41 days" appears to be wrong, unsupported, or from a different source — it does not match this paper's reported 0.96 months (≈29 days).

### Cited Findings
- NSCLC: "median time to NGS was 0.96 (95% CI, 0.93-0.96) months" — the fastest of the five cancer types — [JAMA Network Open, full-text retrieval](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- Metastatic prostate cancer: "median time to NGS was 11.0 (95% CI, 9.8-11.0) months" — the slowest of the five cancer types, matching the claimed "~11.0 months" figure exactly — [JAMA Network Open, full-text retrieval](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- Metastatic breast cancer: 2.9 months (95% CI, 2.7–3.1) — [JAMA Network Open, full-text retrieval](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- Metastatic colorectal cancer: 1.5 months (95% CI, 1.5–1.6) — [JAMA Network Open, full-text retrieval](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- Metastatic pancreatic cancer: 1.3 months (95% CI, 1.2–1.3) — [JAMA Network Open, full-text retrieval](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2847397)
- 0.96 months converts to approximately 29.2 days (0.96 × ~30.44 days/month), not 41 days — arithmetic conversion, no independent source needed

### Inferences
- "41 days" is most likely either a transcription/misremembering error, a figure from a different (uncited) paper or press summary, or a conflation with some other statistic (e.g., a mean rather than median, or a different cohort/subgroup, or possibly a different cancer-testing study entirely). None of the sources retrieved in this research — including the primary JAMA Network Open page itself — state "41 days" for NSCLC anywhere.
- The correct, citable figure for a slide is: **"NSCLC: median time to NGS = 0.96 months (95% CI, 0.93–0.96), i.e., under 1 month"** — using the paper's own units (months), not a day-count conversion, since the paper does not itself report days.
- The directional claim (NSCLC fastest, metastatic prostate slowest at ~11.0 months) is fully supported and should be kept as-is.

### Gaps
- I retrieved these month-level figures via a content-extraction pass over the full JAMA Network Open article page rather than viewing a raw PDF/screenshot of the results table myself; I was not able to additionally cross-check these exact figures against a second independent source (e.g., a PMC mirror or a secondary outlet that reproduces the full results table) in the time available. The figures are internally consistent (NSCLC fastest, prostate slowest, matching the claimed ~11.0 months exactly) and came directly from the publisher's own page, so confidence is high, but a presenter with time should pull a screenshot of the actual results table/figure directly from jamanetwork.com before finalizing slide text.

## What does the actual STAT News article say, and does it match the claimed April 7, 2026 date and "over 40%" quote?

### Takeaway
The STAT News article is real, and its publication date (April 7, 2026) matches the claim exactly. However, its exact headline and visible text do **not** contain an "over 40%" figure — the freely visible (non-paywalled) portion instead says "roughly half of patients in the cohort received genetic sequencing" (implying roughly half did not). The article is paywalled as a STAT+ exclusive beyond the introductory paragraphs, so the full text — which might state "over 40%" somewhere further down — could not be retrieved and verified.

### Cited Findings
- Confirmed real article, found directly on statnews.com: title "Many cancer patients don't get genomic tests to guide treatment, study finds," by Angus Chen (STAT's cancer reporter), published April 7, 2026 — [STAT News](https://www.statnews.com/2026/04/07/genomic-tests-needed-for-targeted-cancer-treatments/)
- The visible (non-paywalled) lede states: "For most cancers in the study, roughly half of patients in the cohort received genetic sequencing" — i.e., the article's own framing in its open portion is "roughly half" tested (so roughly half untested), not an explicit "over 40%" figure — [STAT News](https://www.statnews.com/2026/04/07/genomic-tests-needed-for-targeted-cancer-treatments/)
- The article is marked "exclusive to STAT+ subscribers" and the retrievable content cuts off after the introductory paragraphs, so the full body text (which may include a more precise percentage further down) was not accessible — [STAT News](https://www.statnews.com/2026/04/07/genomic-tests-needed-for-targeted-cancer-treatments/)
- The article covers the same JAMA Network Open study (patients with five metastatic/advanced cancers; disparities by income, insurance type [Medicare/Medicaid], and race/ethnicity [Black, Hispanic] in testing access) — consistent with the JAMA paper verified above — [STAT News](https://www.statnews.com/2026/04/07/genomic-tests-needed-for-targeted-cancer-treatments/)

### Inferences
- The claimed STAT article is real and correctly dated, but the specific "over 40%" wording/quote could not be confirmed from the accessible portion of the text. Given that the JAMA paper's own aggregate (per the derived ~41.3% figure, corroborated by Medscape) is "roughly 41%," it is plausible that STAT's full (paywalled) article states "over 40%" somewhere in the body — this would be consistent with, not contradicted by, what is visible — but it is not something I could verify directly.
- For precision in a citation, it would be safer either to (a) obtain STAT+ access and quote the exact sentence, or (b) cite only what's independently verifiable: the JAMA paper's per-cancer-type and ~41% pooled non-testing figures, attributing the "over 40%" framing to secondary reporting broadly (Medscape, STAT) rather than asserting it is a verbatim STAT quote without having read the full paywalled text.

### Gaps
- Full STAT+ paywalled text was not accessible, so the exact "over 40%" sentence (if present) could not be confirmed verbatim. A presenter with a STAT+ subscription should pull the exact quote directly before using it as a verbatim citation.
