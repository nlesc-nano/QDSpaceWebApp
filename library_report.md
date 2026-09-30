# Library index report

- legacy structure groups: 174
- same start, other functional (merged): 23
- merged with builder twins: 20
- meta.yaml written: 0
- id collisions (suffixed): 48
- errors: 2

| family | material | structures | builder | dft | builder+dft | reconstructed | centres |
|---|---|---|---|---|---|---|---|
| ABX3 | CsPbBr3 | 11 | 5 | 5 | 1 | 0 | Cs:11 |
| ABX3 | CsPbCl3 | 11 | 6 | 5 | 0 | 0 | Cs:10, interstitial:1 |
| ABX3 | CsPbI3 | 8 | 6 | 2 | 0 | 0 | Cs:6, Pb:2 |
| II-VI | CdS | 71 | 66 | 2 | 3 | 33 | Cd:35, S:36 |
| II-VI | CdSe | 79 | 66 | 10 | 3 | 35 | Cd:35, Se:44 |
| II-VI | CdTe | 74 | 69 | 5 | 0 | 33 | Cd:35, Te:39 |
| II-VI | HgS | 75 | 70 | 2 | 3 | 35 | Hg:37, S:38 |
| II-VI | HgSe | 75 | 70 | 2 | 3 | 35 | Hg:37, Se:38 |
| II-VI | HgTe | 78 | 73 | 5 | 0 | 35 | Hg:37, Te:41 |
| II-VI | ZnS | 73 | 68 | 4 | 1 | 33 | S:38, Zn:35 |
| II-VI | ZnSe | 72 | 67 | 3 | 2 | 33 | Se:37, Zn:35 |
| II-VI | ZnTe | 71 | 66 | 2 | 3 | 33 | Te:36, Zn:35 |
| II-VI@II-VI | CdS_ZnS | 1 | 0 | 1 | 0 | 0 | S:1 |
| II-VI@II-VI | CdS_ZnSe | 1 | 0 | 1 | 0 | 0 | S:1 |
| II-VI@II-VI | CdSe_CdS | 1 | 0 | 1 | 0 | 0 | Se:1 |
| II-VI@II-VI | CdSe_ZnS | 1 | 0 | 1 | 0 | 0 | Se:1 |
| II-VI@II-VI | CdSe_ZnSe | 1 | 0 | 1 | 0 | 0 | Se:1 |
| II-VI@II-VI | ZnSe_CdS | 1 | 0 | 1 | 0 | 0 | Se:1 |
| II-VI@II-VI | ZnSe_ZnS | 1 | 0 | 1 | 0 | 0 | Se:1 |
| III-V | GaAs | 40 | 30 | 10 | 0 | 10 | As:16, Ga:24 |
| III-V | GaP | 40 | 30 | 10 | 0 | 10 | Ga:24, P:16 |
| III-V | GaSb | 38 | 30 | 8 | 0 | 10 | Ga:22, Sb:16 |
| III-V | InAs | 45 | 30 | 15 | 0 | 10 | As:16, In:29 |
| III-V | InP | 43 | 29 | 13 | 1 | 10 | In:27, P:16 |
| III-V | InSb | 39 | 30 | 9 | 0 | 10 | In:23, Sb:16 |
| IV-VI | PbS | 31 | 27 | 4 | 0 | 0 | Pb:17, S:14 |
| IV-VI | PbSe | 30 | 27 | 3 | 0 | 0 | Pb:16, Se:14 |
| IV-VI | PbTe | 30 | 27 | 3 | 0 | 0 | Pb:16, Te:14 |

## Merged with builder twins

- ABX3/CsPbBr3 -> CsPbBr3-Cs-Cs20Pb8Br36-clean (Δr̄ = 0.000 Å)
- II-VI/CdS/HLE17/12ang -> CdS-S-Cd16S13Cl6-clean (Δr̄ = 0.346 Å)
- II-VI/CdS/HLE17/20ang -> CdS-S-Cd68S55Cl26-clean (Δr̄ = 0.251 Å)
- II-VI/CdS/HLE17/28ang -> CdS-S-Cd176S147Cl58-clean (Δr̄ = 0.466 Å)
- II-VI/CdSe/HLE17/12ang -> CdSe-Se-Cd16Se13Cl6-clean (Δr̄ = 0.388 Å)
- II-VI/CdSe/HLE17/20ang -> CdSe-Se-Cd68Se55Cl26-clean (Δr̄ = 0.265 Å)
- II-VI/CdSe/HLE17/28ang -> CdSe-Se-Cd176Se147Cl58-clean (Δr̄ = 0.209 Å)
- II-VI/HgS/HLE17/12ang, II-VI/HgS/PBE/12ang -> HgS-S-Hg16S13Cl6-clean (Δr̄ = 0.355 Å)
- II-VI/HgS/HLE17/20ang, II-VI/HgS/PBE/20ang -> HgS-S-Hg68S55Cl26-clean (Δr̄ = 0.230 Å)
- II-VI/HgS/HLE17/28ang, II-VI/HgS/PBE/28ang -> HgS-S-Hg176S147Cl58-clean (Δr̄ = 0.391 Å)
- II-VI/HgSe/HLE17/12ang, II-VI/HgSe/PBE/12ang -> HgSe-Se-Hg16Se13Cl6-clean (Δr̄ = 0.397 Å)
- II-VI/HgSe/HLE17/20ang, II-VI/HgSe/PBE/20ang -> HgSe-Se-Hg68Se55Cl26-clean (Δr̄ = 0.289 Å)
- II-VI/HgSe/HLE17/28ang, II-VI/HgSe/PBE/28ang -> HgSe-Se-Hg176Se147Cl58-clean (Δr̄ = 0.216 Å)
- II-VI/ZnS/HLE17/12ang -> ZnS-S-Zn16S13Cl6-clean (Δr̄ = 0.380 Å)
- II-VI/ZnSe/HLE17/12ang -> ZnSe-Se-Zn16Se13Cl6-clean (Δr̄ = 0.331 Å)
- II-VI/ZnSe/HLE17/20ang -> ZnSe-Se-Zn68Se55Cl26-clean (Δr̄ = 0.423 Å)
- II-VI/ZnTe/HLE17/12ang -> ZnTe-Te-Zn16Te13Cl6-clean (Δr̄ = 0.383 Å)
- II-VI/ZnTe/HLE17/20ang -> ZnTe-Te-Zn68Te55Cl26-clean (Δr̄ = 0.256 Å)
- II-VI/ZnTe/HLE17/28ang -> ZnTe-Te-Zn176Te147Cl58-clean (Δr̄ = 0.215 Å)
- III-V/InP/old/InP_84.xyz -> InP-In-In31P20Cl33-clean (Δr̄ = 0.492 Å)

## Same start geometry, other functional

- II-VI/HgS/PBE/12ang -> II-VI/HgS/HLE17/12ang
- II-VI/HgS/PBE/20ang -> II-VI/HgS/HLE17/20ang
- II-VI/HgS/PBE/28ang -> II-VI/HgS/HLE17/28ang
- II-VI/HgS/PBE/34ang -> II-VI/HgS/HLE17/34ang
- II-VI/HgS/PBE/40ang -> II-VI/HgS/HLE17/40ang
- II-VI/HgSe/PBE/12ang -> II-VI/HgSe/HLE17/12ang
- II-VI/HgSe/PBE/20ang -> II-VI/HgSe/HLE17/20ang
- II-VI/HgSe/PBE/28ang -> II-VI/HgSe/HLE17/28ang
- II-VI/HgSe/PBE/34ang -> II-VI/HgSe/HLE17/34ang
- II-VI/HgSe/PBE/40ang -> II-VI/HgSe/HLE17/40ang
- II-VI/HgTe/PBE/12ang -> II-VI/HgTe/HLE17/12ang
- II-VI/HgTe/PBE/20ang -> II-VI/HgTe/HLE17/20ang
- II-VI/HgTe/PBE/28ang -> II-VI/HgTe/HLE17/28ang
- II-VI/HgTe/PBE/34ang -> II-VI/HgTe/HLE17/34ang
- II-VI/HgTe/PBE/40ang -> II-VI/HgTe/HLE17/40ang
- III-V/InP/old/InP_1116.xyz -> III-V/InP/HLE17/35ang
- III-V/InP/old/InP_1116_ord.xyz -> III-V/InP/HLE17/35ang
- III-V/InP/old/InP_1847.xyz -> III-V/InP/HLE17/41ang
- III-V/InP/old/InP_282.xyz -> III-V/InP/HLE17/23ang
- III-V/InP/old/InP_2839.xyz -> III-V/InP/HLE17/47ang
- III-V/InP/old/InP_608.xyz -> III-V/InP/HLE17/29ang
- III-V/InP/old/InP_84_ord.xyz -> III-V/InP/HLE17/18ang
- III-V/InP/old/mol.xyz -> III-V/InP/HLE17/47ang

## Id collisions

- CsPbBr3-Cs-Cs20Pb8Br36-clean -> CsPbBr3-Cs-Cs20Pb8Br36-clean-v2 (ABX3/CsPbBr3/HLE17/12ang)
- CsPbBr3-Cs-Cs20Pb8Br36-clean -> CsPbBr3-Cs-Cs20Pb8Br36-clean-v3 (ABX3/CsPbBr3/HLE17/16ang)
- CsPbBr3-Cs-Cs112Pb64Br240-clean -> CsPbBr3-Cs-Cs112Pb64Br240-clean-v2 (ABX3/CsPbBr3/HLE17/24ang)
- CsPbBr3-Cs-Cs324Pb216Br756-clean -> CsPbBr3-Cs-Cs324Pb216Br756-clean-v2 (ABX3/CsPbBr3/HLE17/37ang)
- CsPbBr3-Cs-Cs704Pb512Br1728-clean -> CsPbBr3-Cs-Cs704Pb512Br1728-clean-v2 (ABX3/CsPbBr3/HLE17/48ang)
- CsPbCl3-Cs-Cs20Pb8Cl36-clean -> CsPbCl3-Cs-Cs20Pb8Cl36-clean-v2 (ABX3/CsPbCl3/HLE17/12ang)
- CsPbCl3-Cs-Cs112Pb64Cl240-clean -> CsPbCl3-Cs-Cs112Pb64Cl240-clean-v2 (ABX3/CsPbCl3/HLE17/24ang)
- CdSe-Se-Cd176Se147Cl58-clean -> CdSe-Se-Cd176Se147Cl58-clean-v2 (II-VI/CdSe/old/Cd176_OPT.xyz)
- CdSe-Se-Cd176Se147Cl58-clean -> CdSe-Se-Cd176Se147Cl58-clean-v3 (II-VI/CdSe/old/Cd176_pure_OPT.xyz)
- CdSe-Se-Cd68Se55Cl26-clean -> CdSe-Se-Cd68Se55Cl26-clean-v2 (II-VI/CdSe/old/Cd68_OPT.xyz)
- CdSe-Se-Cd360Se309Cl102-clean -> CdSe-Se-Cd360Se309Cl102-clean-v2 (II-VI/CdSe/old/CdSe_771-pos-1.xyz)
- CdSe-Se-Cd310Se281Cl58-clean -> CdSe-Se-Cd310Se281Cl58-clean-v2 (II-VI/CdSe/old/mol.xyz)
- CdTe-Te-Cd16Te13Cl6-clean -> CdTe-Te-Cd16Te13Cl6-clean-v2 (II-VI/CdTe/HLE17/12ang)
- CdTe-Te-Cd68Te55Cl26-clean -> CdTe-Te-Cd68Te55Cl26-clean-v2 (II-VI/CdTe/HLE17/20ang)
- CdTe-Te-Cd176Te147Cl58-clean -> CdTe-Te-Cd176Te147Cl58-clean-v2 (II-VI/CdTe/HLE17/28ang)
- HgTe-Te-Hg16Te13Cl6-clean -> HgTe-Te-Hg16Te13Cl6-clean-v2 (II-VI/HgTe/HLE17/12ang, II-VI/HgTe/PBE/12ang)
- HgTe-Te-Hg68Te55Cl26-clean -> HgTe-Te-Hg68Te55Cl26-clean-v2 (II-VI/HgTe/HLE17/20ang, II-VI/HgTe/PBE/20ang)
- HgTe-Te-Hg176Te147Cl58-clean -> HgTe-Te-Hg176Te147Cl58-clean-v2 (II-VI/HgTe/HLE17/28ang, II-VI/HgTe/PBE/28ang)
- ZnS-S-Zn68S55Cl26-clean -> ZnS-S-Zn68S55Cl26-clean-v2 (II-VI/ZnS/HLE17/20ang)
- ZnS-S-Zn176S147Cl58-clean -> ZnS-S-Zn176S147Cl58-clean-v2 (II-VI/ZnS/HLE17/28ang)
- ZnSe-Se-Zn176Se147Cl58-clean -> ZnSe-Se-Zn176Se147Cl58-clean-v2 (II-VI/ZnSe/HLE17/28ang)
- GaAs-Ga-Ga31As20Cl33-clean -> GaAs-Ga-Ga31As20Cl33-clean-v2 (III-V/GaAs/HLE17/18ang)
- GaAs-Ga-Ga107As73Cl102-clean -> GaAs-Ga-Ga107As73Cl102-clean-v2 (III-V/GaAs/HLE17/23ang)
- GaAs-Ga-Ga249As194Cl165-clean -> GaAs-Ga-Ga249As194Cl165-clean-v2 (III-V/GaAs/HLE17/29ang)
- GaAs-Ga-Ga477As396Cl243-clean -> GaAs-Ga-Ga477As396Cl243-clean-v2 (III-V/GaAs/HLE17/35ang)
- GaP-Ga-Ga31P20Cl33-clean -> GaP-Ga-Ga31P20Cl33-clean-v2 (III-V/GaP/HLE17/18ang)
- GaP-Ga-Ga249P194Cl165-clean -> GaP-Ga-Ga249P194Cl165-clean-v2 (III-V/GaP/HLE17/29ang)
- GaP-Ga-Ga477P396Cl243-clean -> GaP-Ga-Ga477P396Cl243-clean-v2 (III-V/GaP/HLE17/35ang)
- GaP-Ga-Ga790P688Cl348Zn21-clean -> GaP-Ga-Ga790P688Cl348Zn21-clean-v2 (III-V/GaP/HLE17/41ang)
- GaSb-Ga-Ga31Sb20Cl33-clean -> GaSb-Ga-Ga31Sb20Cl33-clean-v2 (III-V/GaSb/HLE17/18ang)
- GaSb-Ga-Ga249Sb194Cl165-clean -> GaSb-Ga-Ga249Sb194Cl165-clean-v2 (III-V/GaSb/HLE17/30ang)
- InAs-In-In31As20Cl33-clean -> InAs-In-In31As20Cl33-clean-v2 (III-V/InAs/HLE17/18ang)
- InAs-In-In107As73Cl102-clean -> InAs-In-In107As73Cl102-clean-v2 (III-V/InAs/HLE17/24ang)
- InAs-In-In249As194Cl165-clean -> InAs-In-In249As194Cl165-clean-v2 (III-V/InAs/HLE17/30ang)
- InAs-In-In477As396Cl243-clean -> InAs-In-In477As396Cl243-clean-v2 (III-V/InAs/HLE17/36ang)
- InAs-In-In790As688Cl348Zn21-clean -> InAs-In-In790As688Cl348Zn21-clean-v2 (III-V/InAs/HLE17/42ang)
- InAs-In-In477As396Cl243-clean -> InAs-In-In477As396Cl243-clean-v3 (III-V/InAs/OLD/InAs1116_GeoOpt-pos-1.xyz)
- InAs-In-In107As73Cl102-clean -> InAs-In-In107As73Cl102-clean-v3 (III-V/InAs/OLD/InAs282_GeoOpt-pos-1.xyz)
- InAs-In-In249As194Cl165-clean -> InAs-In-In249As194Cl165-clean-v3 (III-V/InAs/OLD/InAs608_GeoOpt-pos-1.xyz)
- InP-In-In477P396Cl243-clean -> InP-In-In477P396Cl243-clean-v2 (III-V/InP/old/InP_1116_GeoOpt-pos-1.xyz)
- InP-In-In790P688Cl348Zn21-clean -> InP-In-In790P688Cl348Zn21-clean-v2 (III-V/InP/old/InP_1847_GeoOpt-pos-1.xyz)
- InP-In-In107P73Cl102-clean -> InP-In-In107P73Cl102-clean-v2 (III-V/InP/old/InP_282_GeoOpt-pos-1.xyz)
- InP-In-In1242P1108Cl460Zn29-clean -> InP-In-In1242P1108Cl460Zn29-clean-v2 (III-V/InP/old/InP_2839_GeoOpt-pos-1.xyz)
- InP-In-In249P194Cl165-clean -> InP-In-In249P194Cl165-clean-v2 (III-V/InP/old/InP_608_GeoOpt-pos-1.xyz)
- InSb-In-In31Sb20Cl33-clean -> InSb-In-In31Sb20Cl33-clean-v2 (III-V/InSb/HLE17/18ang)
- InSb-In-In249Sb194Cl165-clean -> InSb-In-In249Sb194Cl165-clean-v2 (III-V/InSb/HLE17/32ang)
- InSb-In-In477Sb396Cl243-clean -> InSb-In-In477Sb396Cl243-clean-v2 (III-V/InSb/HLE17/39ang)
- InSb-In-In790Sb688Cl348Zn21-clean -> InSb-In-In790Sb688Cl348Zn21-clean-v2 (III-V/InSb/HLE17/45ang)

## Errors

- III-V/InAs/OLD/InAs1847_GeoOpt-pos-1.xyz: IndexError: list index out of range
- III-V/InAs/OLD/InAs2839_GeoOpt-pos-1.xyz: IndexError: list index out of range
