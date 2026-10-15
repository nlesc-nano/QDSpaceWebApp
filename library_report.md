# Library index report

- builder variant-recipe duplicates (merged): 98
- legacy structure groups: 174
- same start, other functional (merged): 23
- merged with builder twins: 38
- meta.yaml written: 0
- id collisions (suffixed): 32
- errors: 2

| family | material | structures | builder | dft | builder+dft | reconstructed | centres |
|---|---|---|---|---|---|---|---|
| ABX3 | CsPbBr3 | 11 | 7 | 0 | 4 | 0 | Cs:6, Pb:5 |
| ABX3 | CsPbCl3 | 14 | 9 | 3 | 2 | 0 | Cs:8, Pb:5, interstitial:1 |
| ABX3 | CsPbI3 | 11 | 9 | 0 | 2 | 0 | Cs:6, Pb:5 |
| II-VI | CdS | 102 | 97 | 2 | 3 | 64 | Cd:51, S:51 |
| II-VI | CdSe | 278 | 268 | 6 | 4 | 180 | Cd:139, Se:139 |
| II-VI | CdTe | 105 | 100 | 5 | 0 | 64 | Cd:51, Te:54 |
| II-VI | HgS | 109 | 104 | 2 | 3 | 69 | Hg:54, S:55 |
| II-VI | HgSe | 109 | 104 | 2 | 3 | 69 | Hg:54, Se:55 |
| II-VI | HgTe | 112 | 107 | 5 | 0 | 69 | Hg:54, Te:58 |
| II-VI | ZnS | 104 | 99 | 4 | 1 | 64 | S:53, Zn:51 |
| II-VI | ZnSe | 103 | 98 | 3 | 2 | 64 | Se:52, Zn:51 |
| II-VI | ZnTe | 103 | 98 | 2 | 3 | 65 | Te:51, Zn:52 |
| II-VI@II-VI | CdS_ZnS | 1 | 0 | 1 | 0 | 0 | S:1 |
| II-VI@II-VI | CdS_ZnSe | 1 | 0 | 1 | 0 | 0 | S:1 |
| II-VI@II-VI | CdSe_CdS | 1 | 0 | 1 | 0 | 0 | Se:1 |
| II-VI@II-VI | CdSe_ZnS | 1 | 0 | 1 | 0 | 0 | Se:1 |
| II-VI@II-VI | CdSe_ZnSe | 1 | 0 | 1 | 0 | 0 | Se:1 |
| II-VI@II-VI | ZnSe_CdS | 1 | 0 | 1 | 0 | 0 | Se:1 |
| II-VI@II-VI | ZnSe_ZnS | 1 | 0 | 1 | 0 | 0 | Se:1 |
| III-V | AlAs | 52 | 52 | 0 | 0 | 22 | Al:27, As:25 |
| III-V | AlP | 52 | 52 | 0 | 0 | 22 | Al:27, P:25 |
| III-V | AlSb | 52 | 52 | 0 | 0 | 22 | Al:27, Sb:25 |
| III-V | GaAs | 61 | 51 | 9 | 1 | 22 | As:26, Ga:35 |
| III-V | GaP | 61 | 51 | 9 | 1 | 22 | Ga:35, P:26 |
| III-V | GaSb | 59 | 51 | 7 | 1 | 22 | Ga:33, Sb:26 |
| III-V | InAs | 66 | 51 | 14 | 1 | 22 | As:26, In:40 |
| III-V | InP | 65 | 51 | 13 | 1 | 22 | In:39, P:26 |
| III-V | InSb | 60 | 51 | 8 | 1 | 22 | In:34, Sb:26 |
| IV-VI | PbS | 29 | 25 | 4 | 0 | 0 | Pb:16, S:13 |
| IV-VI | PbSe | 28 | 25 | 3 | 0 | 0 | Pb:15, Se:13 |
| IV-VI | PbTe | 28 | 25 | 3 | 0 | 0 | Pb:15, Te:13 |

## Builder variant-recipe duplicates

- cds_zb_swap.yaml: CdS-Cd-Cd55S28Cl54-reconstructed -> CdS-Cd-Cd55S28Cl54-reconstructed
- cds_zb_swap.yaml: CdS-S-Cd92S55Cl74-reconstructed -> CdS-S-Cd92S55Cl74-reconstructed
- cdse_zb_swap.yaml: CdSe-Cd-Cd55Se28Cl54-reconstructed -> CdSe-Cd-Cd55Se28Cl54-reconstructed
- cdse_zb_swap.yaml: CdSe-Se-Cd92Se55Cl74-reconstructed -> CdSe-Se-Cd92Se55Cl74-reconstructed
- cdse_wz_swap.yaml: CdSe-Cd-Cd16Se7Cl18-reconstructed -> CdSe-Cd-Cd16Se7Cl18-reconstructed
- cdse_wz_swap.yaml: CdSe-Cd-Cd21Se13Cl16-reconstructed -> CdSe-Cd-Cd21Se13Cl16-reconstructed
- cdte_zb_swap.yaml: CdTe-Cd-Cd55Te28Cl54-reconstructed -> CdTe-Cd-Cd55Te28Cl54-reconstructed
- cdte_zb_swap.yaml: CdTe-Te-Cd92Te55Cl74-reconstructed -> CdTe-Te-Cd92Te55Cl74-reconstructed
- hgs_zb_swap.yaml: HgS-Hg-Hg55S28Cl54-reconstructed -> HgS-Hg-Hg55S28Cl54-reconstructed
- hgse_zb_swap.yaml: HgSe-Hg-Hg55Se28Cl54-reconstructed -> HgSe-Hg-Hg55Se28Cl54-reconstructed
- hgte_zb_swap.yaml: HgTe-Hg-Hg55Te28Cl54-reconstructed -> HgTe-Hg-Hg55Te28Cl54-reconstructed
- zns_zb_swap.yaml: ZnS-S-Zn92S55Cl74-reconstructed -> ZnS-S-Zn92S55Cl74-reconstructed
- zns_zb_swap.yaml: ZnS-Zn-Zn55S28Cl54-reconstructed -> ZnS-Zn-Zn55S28Cl54-reconstructed
- znse_zb_swap.yaml: ZnSe-Se-Zn92Se55Cl74-reconstructed -> ZnSe-Se-Zn92Se55Cl74-reconstructed
- znse_zb_swap.yaml: ZnSe-Zn-Zn55Se28Cl54-reconstructed -> ZnSe-Zn-Zn55Se28Cl54-reconstructed
- znte_zb_swap.yaml: ZnTe-Te-Zn92Te55Cl74-reconstructed -> ZnTe-Te-Zn92Te55Cl74-reconstructed
- znte_zb_swap.yaml: ZnTe-Zn-Zn55Te28Cl54-reconstructed -> ZnTe-Zn-Zn55Te28Cl54-reconstructed
- alas_zb_100.yaml: AlAs-Al-Al161As139Cl66-clean -> AlAs-Al-Al161As139Cl66-clean
- alas_zb_100.yaml: AlAs-Al-Al165As139Cl78-clean -> AlAs-Al-Al165As139Cl78-clean
- alas_zb_100.yaml: AlAs-Al-Al35As26Cl27-clean -> AlAs-Al-Al35As26Cl27-clean
- alas_zb_100.yaml: AlAs-Al-Al439As398Cl123-clean -> AlAs-Al-Al439As398Cl123-clean
- alas_zb_100.yaml: AlAs-As-Al348As314Cl102-clean -> AlAs-As-Al348As314Cl102-clean
- alas_zb_100.yaml: AlAs-As-Al360As320Cl120-clean -> AlAs-As-Al360As320Cl120-clean
- alas_zb_swap.yaml: AlAs-Al-Al31As16Cl45-reconstructed -> AlAs-Al-Al31As16Cl45-reconstructed
- alas_zb_swap.yaml: AlAs-As-Al40As23Cl51-reconstructed -> AlAs-As-Al40As23Cl51-reconstructed
- alas_zb_swap.yaml: AlAs-As-Al52As31Cl63-reconstructed -> AlAs-As-Al52As31Cl63-reconstructed
- alp_zb_100.yaml: AlP-Al-Al161P139Cl66-clean -> AlP-Al-Al161P139Cl66-clean
- alp_zb_100.yaml: AlP-Al-Al165P139Cl78-clean -> AlP-Al-Al165P139Cl78-clean
- alp_zb_100.yaml: AlP-Al-Al35P26Cl27-clean -> AlP-Al-Al35P26Cl27-clean
- alp_zb_100.yaml: AlP-Al-Al439P398Cl123-clean -> AlP-Al-Al439P398Cl123-clean
- alp_zb_100.yaml: AlP-P-Al348P314Cl102-clean -> AlP-P-Al348P314Cl102-clean
- alp_zb_100.yaml: AlP-P-Al360P320Cl120-clean -> AlP-P-Al360P320Cl120-clean
- alp_zb_swap.yaml: AlP-Al-Al31P16Cl45-reconstructed -> AlP-Al-Al31P16Cl45-reconstructed
- alp_zb_swap.yaml: AlP-P-Al40P23Cl51-reconstructed -> AlP-P-Al40P23Cl51-reconstructed
- alp_zb_swap.yaml: AlP-P-Al52P31Cl63-reconstructed -> AlP-P-Al52P31Cl63-reconstructed
- alsb_zb_100.yaml: AlSb-Al-Al161Sb139Cl66-clean -> AlSb-Al-Al161Sb139Cl66-clean
- alsb_zb_100.yaml: AlSb-Al-Al165Sb139Cl78-clean -> AlSb-Al-Al165Sb139Cl78-clean
- alsb_zb_100.yaml: AlSb-Al-Al35Sb26Cl27-clean -> AlSb-Al-Al35Sb26Cl27-clean
- alsb_zb_100.yaml: AlSb-Al-Al439Sb398Cl123-clean -> AlSb-Al-Al439Sb398Cl123-clean
- alsb_zb_100.yaml: AlSb-Sb-Al348Sb314Cl102-clean -> AlSb-Sb-Al348Sb314Cl102-clean
- alsb_zb_100.yaml: AlSb-Sb-Al360Sb320Cl120-clean -> AlSb-Sb-Al360Sb320Cl120-clean
- alsb_zb_swap.yaml: AlSb-Al-Al31Sb16Cl45-reconstructed -> AlSb-Al-Al31Sb16Cl45-reconstructed
- alsb_zb_swap.yaml: AlSb-Sb-Al40Sb23Cl51-reconstructed -> AlSb-Sb-Al40Sb23Cl51-reconstructed
- alsb_zb_swap.yaml: AlSb-Sb-Al52Sb31Cl63-reconstructed -> AlSb-Sb-Al52Sb31Cl63-reconstructed
- gaas_zb_100.yaml: GaAs-As-Ga348As314Cl102-clean -> GaAs-As-Ga348As314Cl102-clean
- gaas_zb_100.yaml: GaAs-As-Ga360As320Cl120-clean -> GaAs-As-Ga360As320Cl120-clean
- gaas_zb_100.yaml: GaAs-Ga-Ga161As139Cl66-clean -> GaAs-Ga-Ga161As139Cl66-clean
- gaas_zb_100.yaml: GaAs-Ga-Ga165As139Cl78-clean -> GaAs-Ga-Ga165As139Cl78-clean
- gaas_zb_100.yaml: GaAs-Ga-Ga35As26Cl27-clean -> GaAs-Ga-Ga35As26Cl27-clean
- gaas_zb_100.yaml: GaAs-Ga-Ga439As398Cl123-clean -> GaAs-Ga-Ga439As398Cl123-clean
- gaas_zb_swap.yaml: GaAs-As-Ga40As23Cl51-reconstructed -> GaAs-As-Ga40As23Cl51-reconstructed
- gaas_zb_swap.yaml: GaAs-As-Ga52As31Cl63-reconstructed -> GaAs-As-Ga52As31Cl63-reconstructed
- gaas_zb_swap.yaml: GaAs-Ga-Ga31As16Cl45-reconstructed -> GaAs-Ga-Ga31As16Cl45-reconstructed
- gap_zb_100.yaml: GaP-Ga-Ga161P139Cl66-clean -> GaP-Ga-Ga161P139Cl66-clean
- gap_zb_100.yaml: GaP-Ga-Ga165P139Cl78-clean -> GaP-Ga-Ga165P139Cl78-clean
- gap_zb_100.yaml: GaP-Ga-Ga35P26Cl27-clean -> GaP-Ga-Ga35P26Cl27-clean
- gap_zb_100.yaml: GaP-Ga-Ga439P398Cl123-clean -> GaP-Ga-Ga439P398Cl123-clean
- gap_zb_100.yaml: GaP-P-Ga348P314Cl102-clean -> GaP-P-Ga348P314Cl102-clean
- gap_zb_100.yaml: GaP-P-Ga360P320Cl120-clean -> GaP-P-Ga360P320Cl120-clean
- gap_zb_swap.yaml: GaP-Ga-Ga31P16Cl45-reconstructed -> GaP-Ga-Ga31P16Cl45-reconstructed
- gap_zb_swap.yaml: GaP-P-Ga40P23Cl51-reconstructed -> GaP-P-Ga40P23Cl51-reconstructed
- gap_zb_swap.yaml: GaP-P-Ga52P31Cl63-reconstructed -> GaP-P-Ga52P31Cl63-reconstructed
- gasb_zb_100.yaml: GaSb-Ga-Ga161Sb139Cl66-clean -> GaSb-Ga-Ga161Sb139Cl66-clean
- gasb_zb_100.yaml: GaSb-Ga-Ga165Sb139Cl78-clean -> GaSb-Ga-Ga165Sb139Cl78-clean
- gasb_zb_100.yaml: GaSb-Ga-Ga35Sb26Cl27-clean -> GaSb-Ga-Ga35Sb26Cl27-clean
- gasb_zb_100.yaml: GaSb-Ga-Ga439Sb398Cl123-clean -> GaSb-Ga-Ga439Sb398Cl123-clean
- gasb_zb_100.yaml: GaSb-Sb-Ga348Sb314Cl102-clean -> GaSb-Sb-Ga348Sb314Cl102-clean
- gasb_zb_100.yaml: GaSb-Sb-Ga360Sb320Cl120-clean -> GaSb-Sb-Ga360Sb320Cl120-clean
- gasb_zb_swap.yaml: GaSb-Ga-Ga31Sb16Cl45-reconstructed -> GaSb-Ga-Ga31Sb16Cl45-reconstructed
- gasb_zb_swap.yaml: GaSb-Sb-Ga40Sb23Cl51-reconstructed -> GaSb-Sb-Ga40Sb23Cl51-reconstructed
- gasb_zb_swap.yaml: GaSb-Sb-Ga52Sb31Cl63-reconstructed -> GaSb-Sb-Ga52Sb31Cl63-reconstructed
- inas_zb_100.yaml: InAs-As-In348As314Cl102-clean -> InAs-As-In348As314Cl102-clean
- inas_zb_100.yaml: InAs-As-In360As320Cl120-clean -> InAs-As-In360As320Cl120-clean
- inas_zb_100.yaml: InAs-In-In161As139Cl66-clean -> InAs-In-In161As139Cl66-clean
- inas_zb_100.yaml: InAs-In-In165As139Cl78-clean -> InAs-In-In165As139Cl78-clean
- inas_zb_100.yaml: InAs-In-In35As26Cl27-clean -> InAs-In-In35As26Cl27-clean
- inas_zb_100.yaml: InAs-In-In439As398Cl123-clean -> InAs-In-In439As398Cl123-clean
- inas_zb_swap.yaml: InAs-As-In40As23Cl51-reconstructed -> InAs-As-In40As23Cl51-reconstructed
- inas_zb_swap.yaml: InAs-As-In52As31Cl63-reconstructed -> InAs-As-In52As31Cl63-reconstructed
- inas_zb_swap.yaml: InAs-In-In31As16Cl45-reconstructed -> InAs-In-In31As16Cl45-reconstructed
- inp_zb_100.yaml: InP-In-In161P139Cl66-clean -> InP-In-In161P139Cl66-clean
- inp_zb_100.yaml: InP-In-In165P139Cl78-clean -> InP-In-In165P139Cl78-clean
- inp_zb_100.yaml: InP-In-In35P26Cl27-clean -> InP-In-In35P26Cl27-clean
- inp_zb_100.yaml: InP-In-In439P398Cl123-clean -> InP-In-In439P398Cl123-clean
- inp_zb_100.yaml: InP-P-In348P314Cl102-clean -> InP-P-In348P314Cl102-clean
- inp_zb_100.yaml: InP-P-In360P320Cl120-clean -> InP-P-In360P320Cl120-clean
- inp_zb_swap.yaml: InP-In-In31P16Cl45-reconstructed -> InP-In-In31P16Cl45-reconstructed
- inp_zb_swap.yaml: InP-P-In40P23Cl51-reconstructed -> InP-P-In40P23Cl51-reconstructed
- inp_zb_swap.yaml: InP-P-In52P31Cl63-reconstructed -> InP-P-In52P31Cl63-reconstructed
- insb_zb_100.yaml: InSb-In-In161Sb139Cl66-clean -> InSb-In-In161Sb139Cl66-clean
- insb_zb_100.yaml: InSb-In-In165Sb139Cl78-clean -> InSb-In-In165Sb139Cl78-clean
- insb_zb_100.yaml: InSb-In-In35Sb26Cl27-clean -> InSb-In-In35Sb26Cl27-clean
- insb_zb_100.yaml: InSb-In-In439Sb398Cl123-clean -> InSb-In-In439Sb398Cl123-clean
- insb_zb_100.yaml: InSb-Sb-In348Sb314Cl102-clean -> InSb-Sb-In348Sb314Cl102-clean
- insb_zb_100.yaml: InSb-Sb-In360Sb320Cl120-clean -> InSb-Sb-In360Sb320Cl120-clean
- insb_zb_swap.yaml: InSb-In-In31Sb16Cl45-reconstructed -> InSb-In-In31Sb16Cl45-reconstructed
- insb_zb_swap.yaml: InSb-Sb-In40Sb23Cl51-reconstructed -> InSb-Sb-In40Sb23Cl51-reconstructed
- insb_zb_swap.yaml: InSb-Sb-In52Sb31Cl63-reconstructed -> InSb-Sb-In52Sb31Cl63-reconstructed

## Merged with builder twins

- ABX3/CsPbBr3 -> CsPbBr3-Cs-Cs20Pb8Br36-clean (Δr̄ = 0.000 Å)
- ABX3/CsPbBr3/HLE17/12ang -> CsPbBr3-Cs-Cs20Pb8Br36-clean (Δr̄ = 0.360 Å)
- ABX3/CsPbBr3/HLE17/16ang -> CsPbBr3-Cs-Cs20Pb8Br36-clean (Δr̄ = 0.360 Å)
- ABX3/CsPbBr3/HLE17/24ang -> CsPbBr3-Cs-Cs112Pb64Br240-clean (Δr̄ = 0.314 Å)
- ABX3/CsPbBr3/HLE17/37ang -> CsPbBr3-Cs-Cs324Pb216Br756-clean (Δr̄ = 0.328 Å)
- ABX3/CsPbBr3/HLE17/48ang -> CsPbBr3-Cs-Cs704Pb512Br1728-clean (Δr̄ = 0.417 Å)
- ABX3/CsPbCl3/HLE17/12ang -> CsPbCl3-Cs-Cs20Pb8Cl36-clean (Δr̄ = 0.333 Å)
- ABX3/CsPbCl3/HLE17/24ang -> CsPbCl3-Cs-Cs112Pb64Cl240-clean (Δr̄ = 0.286 Å)
- ABX3/CsPbI3/HLE17/24ang -> CsPbI3-Pb-Cs54Pb27I108-clean (Δr̄ = 0.462 Å)
- ABX3/CsPbI3/HLE17/32ang -> CsPbI3-Pb-Cs200Pb125I450-clean (Δr̄ = 0.662 Å)
- II-VI/CdS/HLE17/12ang -> CdS-S-Cd16S13Cl6-clean (Δr̄ = 0.346 Å)
- II-VI/CdS/HLE17/20ang -> CdS-S-Cd68S55Cl26-clean (Δr̄ = 0.251 Å)
- II-VI/CdS/HLE17/28ang -> CdS-S-Cd176S147Cl58-clean (Δr̄ = 0.466 Å)
- II-VI/CdSe/HLE17/12ang -> CdSe-Se-Cd16Se13Cl6-clean (Δr̄ = 0.388 Å)
- II-VI/CdSe/HLE17/20ang -> CdSe-Se-Cd68Se55Cl26-clean (Δr̄ = 0.265 Å)
- II-VI/CdSe/HLE17/28ang -> CdSe-Se-Cd176Se147Cl58-clean (Δr̄ = 0.209 Å)
- II-VI/CdSe/old/Cd176_OPT.xyz -> CdSe-Se-Cd176Se147Cl58-clean (Δr̄ = 0.209 Å)
- II-VI/CdSe/old/Cd176_pure_OPT.xyz -> CdSe-Se-Cd176Se147Cl58-clean (Δr̄ = 0.209 Å)
- II-VI/CdSe/old/Cd68_OPT.xyz -> CdSe-Se-Cd68Se55Cl26-clean (Δr̄ = 0.261 Å)
- II-VI/CdSe/old/CdSe_771-pos-1.xyz -> CdSe-Se-Cd360Se309Cl102-clean (Δr̄ = 0.176 Å)
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
- III-V/GaAs/HLE17/18ang -> GaAs-Ga-Ga31As20Cl33-clean (Δr̄ = 0.523 Å)
- III-V/GaP/HLE17/18ang -> GaP-Ga-Ga31P20Cl33-clean (Δr̄ = 0.397 Å)
- III-V/GaSb/HLE17/18ang -> GaSb-Ga-Ga31Sb20Cl33-clean (Δr̄ = 0.588 Å)
- III-V/InAs/HLE17/18ang -> InAs-In-In31As20Cl33-clean (Δr̄ = 0.488 Å)
- III-V/InP/old/InP_84.xyz -> InP-In-In31P20Cl33-clean (Δr̄ = 0.492 Å)
- III-V/InSb/HLE17/18ang -> InSb-In-In31Sb20Cl33-clean (Δr̄ = 0.631 Å)

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
- GaAs-Ga-Ga107As73Cl102-clean -> GaAs-Ga-Ga107As73Cl102-clean-v2 (III-V/GaAs/HLE17/23ang)
- GaAs-Ga-Ga249As194Cl165-clean -> GaAs-Ga-Ga249As194Cl165-clean-v2 (III-V/GaAs/HLE17/29ang)
- GaAs-Ga-Ga477As396Cl243-clean -> GaAs-Ga-Ga477As396Cl243-clean-v2 (III-V/GaAs/HLE17/35ang)
- GaP-Ga-Ga249P194Cl165-clean -> GaP-Ga-Ga249P194Cl165-clean-v2 (III-V/GaP/HLE17/29ang)
- GaP-Ga-Ga477P396Cl243-clean -> GaP-Ga-Ga477P396Cl243-clean-v2 (III-V/GaP/HLE17/35ang)
- GaP-Ga-Ga790P688Cl348Zn21-clean -> GaP-Ga-Ga790P688Cl348Zn21-clean-v2 (III-V/GaP/HLE17/41ang)
- GaSb-Ga-Ga249Sb194Cl165-clean -> GaSb-Ga-Ga249Sb194Cl165-clean-v2 (III-V/GaSb/HLE17/30ang)
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
- InSb-In-In249Sb194Cl165-clean -> InSb-In-In249Sb194Cl165-clean-v2 (III-V/InSb/HLE17/32ang)
- InSb-In-In477Sb396Cl243-clean -> InSb-In-In477Sb396Cl243-clean-v2 (III-V/InSb/HLE17/39ang)
- InSb-In-In790Sb688Cl348Zn21-clean -> InSb-In-In790Sb688Cl348Zn21-clean-v2 (III-V/InSb/HLE17/45ang)

## Errors

- III-V/InAs/OLD/InAs1847_GeoOpt-pos-1.xyz: IndexError: list index out of range
- III-V/InAs/OLD/InAs2839_GeoOpt-pos-1.xyz: IndexError: list index out of range
