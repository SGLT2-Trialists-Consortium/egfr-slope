## Column Descriptions

### `synthetic_trial_baseline.csv`

| Variable  | Type      | Description                                                                                           |
| --------- | --------- | ----------------------------------------------------------------------------------------------------- |
| `usubjid` | Character | Unique subject identifier                                                                             |
| `trt01pn` | Integer   | Numeric indicator for treatment assignment; `0` == "Placebo", `1` == "SGLT2i"                         |
| `randfl`  | Character | Randomisation flag (indicator for whether a participant has been randomised)                          |
| `ittfl`   | Character | Intention to treat flag (indicator for whether a participant is in the intention to treat population) |
| `trtfl`   | Character | On-treatment flag (indicator for whether a participant in in the on-treatment population)             |
| `blgfr`   | Integer   | Baseline eGFR value                                                                                   |
| `blglp1`  | Character | Baseline GLP1-RA use                                                                                  |
| `strata`  | Character | Study strata (not required)                                                                           |

------------------------------------------------------------------------------

### `synthetic_trial_follow_up_egfr.csv`

| Variable  | Type      | Description                                                                                           |
| --------- | --------- | ----------------------------------------------------------------------------------------------------- |
| `usubjid` | Character | Unique subject identifier                                                                             |
| `randfl`  | Character | Randomisation flag (indicator for whether a participant has been randomised)                          |
| `ittfl`   | Character | Intention to treat flag (indicator for whether a participant is in the intention to treat population) |
| `trtfl`   | Character | On-treatment flag (indicator for whether a participant in in the on-treatment population)             |
| `anl01fl` | Character | Analysis flag (indicator for whether a value should be used for analysis)                             |
| `aval`    | Double    | Estimated glomerular filtration rate (eGFR) value                                                     |
| `base`    | Integer   | Baseline eGFR value                                                                                   |
| `paramcd` | Character | Parameter code; all `GFRBSCRT`                                                                        |
| `param`   | Character | Parameter value; all `GFR from Creatinine Adjusted for BSA (mL/min/1.73m2)`                           |
| `ady`     | Double    | Days from baseline                                                                                    |
| `avisitn` | Double    | Study visit (weeks from baseline; numeric)                                                            |
| `avisit`  | Character | Study visit (weeks from baseline; character)                                                          |