# QUEST-South — Regional Inversion Application

## Purpose

Regional-scale litho-constrained inversion for the QUEST-South geophysical survey
area, British Columbia. Gravity and magnetic data with geological constraints.

## Status

Data preparation phase. See `regional_QUEST-South/` in the original codebase for
source data. The SimPEG inversion pipeline follows the same pattern as Forrestania.

## Key Files (to be migrated from original)

| File | Purpose |
|------|---------|
| `groom_data.py` | Data preparation from Geoscience BC deliverables |
| `run_joint_inversion.py` | SimPEG joint inversion |
| `*.ini` | C++ inversion configs |

## Dependencies

Same as Forrestania: SimPEG, discretize, cluster_api, lithoseed, all litho-* modules.
