# feature/btj_cloudphys

MPAS-Model v8.4.0 plus the cloud-physics changes used by the JCSDA Global MPAS-JEDI runs on S4
(the "icbl2sv" build), for Skylab and other JEDI users. Each change is a namelist option in
`&physics`, and the defaults reproduce stock MPAS v8.4.0.

| Option | Type | Default | Global runs | What it does |
|---|---|---|---|---|
| `config_icloud_bl` | integer | 0 | 1 | Couples the MYNN subgrid-scale clouds (`cldfrac_bl`, `qc_bl`, `qi_bl`) with RRTMG radiation, as `icloud_bl` in WRF: radiation sees the larger of the resolved and subgrid cloud fraction and condensate. The `cldfrac` written to the output stays the resolved (Xu-Randall) value. Needs `bl_mynn`. |
| `config_gf_shallow_tcap` | real (K) | 0.0 | 273.15 | Grell-Freitas shallow convection: the undilute residual plume stops at the first level above kbcon+tun_lim that is colder than this temperature. 0 = no bound. |
| `config_gf_deep_rhmid` | real | 0.0 | 0.50 | Grell-Freitas deep convection: no deep convection where the pressure-weighted q/qes over 700-400 hPa is below this fraction. 0 = no condition. |
| `config_thompson_snow_vfac` | real | 1.0 | 1.5 | Thompson microphysics: factor on the snow sedimentation (fall) speed at levels with T <= T_0+0.1 K, with the snow sub-step count recomputed for the CFL limit. Deposition and collection terms are unchanged, so only the residence time of snow changes. 1 = original code path. |

Global's namelist also uses `config_physics_suite = 'convection_permitting'`,
`config_convection_scheme = 'suite'`, `config_dt = 60` and 30-minute radiation
(`config_radtlw_interval = config_radtsw_interval = '00:30:00'`), with the standard MPAS-Data tables.

## History

One commit per step of the S4 evaluation, then the namelist options:
icbl (MYNN cloud-radiation coupling) -> icbl2 (coupling moved into the radiation drivers, output
cldfrac stays resolved) -> icbl2s (GF shallow plume-top bound) -> icbl2st (GF deep moisture
condition) -> namelist options in place of a compile-time parameter and environment variables
-> icbl2sv (faster Thompson snow sedimentation at cold levels, `config_thompson_snow_vfac`).

## Reproducing Global's forecasts bit for bit

Checked on a 1-hour forecast on Global's 598,249-cell mesh from 2026-09-07 00Z (128 ranks):
- with Global's option values the forecast is bitwise identical (72 of 72 output variables) to the
  icbl2sv build run with `THOMPSON_SNOW_VFAC=1.5`;
- with Global's values for the first three options and `config_thompson_snow_vfac` at its default,
  it is bitwise identical to the icbl2st build;
- with the defaults it is bitwise identical to unmodified MPAS v8.4.0.

Bitwise agreement needs the same compiler and flags: Intel oneAPI 2024.2.1 (ifort, icx) from
spack-stack 1.9.1, double precision, OpenMP, and Fortran release flags
`-O3 -DNDEBUG -unroll -inline -heap-arrays 32` (the JEDI bundle's Release flags). The same code
built with plain `-O3` gives a different forecast (43 of 72 variables differ after one hour).
