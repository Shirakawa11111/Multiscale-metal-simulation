# evolving_forest_mfp design notes (from design workflow)

## Ledger spec
PER-RECORDED-STEP JSON SCHEMA (one object per REC steps, in result["rows"]; densities in 1/m^2, gamma_m dimensionless):
  istep                 : ExaDiS step index
  gamma_m               : cumulative PRIMARY plastic shear on system m (mobile-only swept area, see below)
  rho_mobile            : length of FAM_m-Burgers arms on plane n_m, not junction-anchored, /vol/b^2
  rho_forest            : length of FAM_f-Burgers arms on plane n_f, not anchored (the live forest)
  rho_junction          : length of b3=b_m+/-b_f-Burgers segments (glissile/Hirth/Lomer product); ==0 for collinear (b3=0)
  rho_residual          : length of FAM_m/FAM_f arms that ARE junction-anchored (degree>=3) or immobilized (the stored residual, incl. collinear residual)
  rho_offsys            : unclassified / ambiguous length (diagnostic; must stay small)
  rho_total             : total network length /vol/b^2
  rho_stored            : rho_junction + rho_residual   (RETAINED immobile product length; THIS is the storage observable)
  rho_removed_events    : cumulative ANNIHILATED length from per-step NEGATIVE rho_total jumps that coincide with a segment-count DROP (i.e. a collision/topology deletion event). This is the direct, mechanism-anchored removal measure.
  rho_removed_balance   : independent cross-check = (rho_built + cum_created) - rho_total - rho_junction, where cum_created is the running sum of POSITIVE rho_total jumps NOT coincident with a seg-count drop (source bow-out + FR multiplication + remesh growth). Demoted to a diagnostic; rho_removed_events is primary.
  mobile_survival       : rho_mobile / peak(rho_mobile)  (carrier-survival gate for the fit window)
  ambiguous_frac        : rho_offsys / rho_total  (ledger-trust gate; reject step if large)
  n_seg                 : segment count (for the seg-drop event detector)

CLASSIFICATION RULE (coarse, no per-segment lineage). For each current segment with endpoints (a,c):
  L = |closest_image(pos[a],pos[c])|  -- ALWAYS PBC-folded (closure/wrap segments otherwise add ~LBOX each).
  PLANE is the PRIMARY discriminator (collinear shares the Burgers vector b_f=+/-b_m, so Burgers cannot
  separate mobile from forest). Burgers-family is the discriminator for NON-collinear pairs; plane corroborates.
  - junction : Burgers || (b_m+/-b_f) and that vector != 0  (collinear has none).
  - collinear case (|cos(b_m,b_f)|>0.94):
        anchored(deg>=3) -> residual ; on-plane n_m -> mobile ; on-plane n_f -> forest ;
        screw/null-plane/off-plane -> velocity tiebreak (|v|>VMOB_THR -> mobile else residual).
  - non-collinear: Burgers==FAM_m -> mobile (if on n_m & not anchored) else residual ;
                   Burgers==FAM_f -> forest (if on n_f & not anchored) else residual ; else offsys.
  Screw segments (plane norm ~ 0) are NEVER normalized (no divide-by-zero); they fall to the velocity tiebreak.

REMOVED / ANNIHILATED LENGTH (two independent definitions, both logged; events is primary):
  (1) rho_removed_events: integrate -d(rho_total) over records where rho_total DROPS AND n_seg drops.
      Rationale: the only sink that removes total length at a topology/collision step is annihilation
      (collinear b_m+b_f->0) or junction zip-into-anchored-residual; remesh deletion of sub-minseg stubs
      also drops n_seg but its magnitude is small and is subtracted out by the same-sense control.
  (2) rho_removed_balance: against an explicit source budget. L_built = total length at the first record;
      cum_created = running sum of positive rho_total jumps with NO seg-drop (sources/multiplication/remesh
      growth). removed = max(0, rho_built + cum_created - rho_total - rho_junction). Used only to bound (1).
  The headline contrast uses rho_removed_events DIFFERENCED against the collinear-SAME control at matched
  gamma_m (post-step): R_remove = d/dgamma[ rho_removed_events(opp) - rho_removed_events(same) ], which
  subtracts the bow-out/remesh/multiplication baseline that is identical between senses.

GAMMA_M (system-m plastic shear), PRIMARY = MOBILE-ONLY SWEPT AREA, accumulated EVERY step:
  positions cached at step_begin (keyed by stable node 'tags', remesh/topology-robust); at step_end, for each
  segment whose Burgers||b_m AND plane||n_m, swept-area kernel n=0.5*cross(r2-r3,r1-r4) (r3,r4 = old endpoints),
  all differences PBC-folded; dEp_m = sym(outer(n,b_m))/V; gamma_m += dEp_m : A2_m (Voigt, factor-2 on shears).
  This is NOT global get_plastic_strain/dEp (which sums ALL slip systems and is fatal for collinear because
  the forest shares b_m). Global dEp resolved on m is accumulated separately as gamma_m_global, used ONLY as a
  loose upper-bound consistency check (gamma_global_consistency = gamma_m/gamma_m_global; flag if << 1).

DENSITY UNITS: vol = |det(cell.h)| in b^3; rho = L_b / vol / b^2 = 1/m^2. L_mf = 1/(b * S_tot) in meters
  (S_tot = S_store + S_remove, units 1/m^2); Lambda = L_mf * sqrt(rho_f_fit), dimensionless.

## Run plan
ONE OS PROCESS = ONE NETWORK = ONE (relax+load) cell  [pyexadis double-free is documented; never loop pairs/controls in one process]. Dispatch as independent HPC jobs, many jobs x few threads (MEMORY: HPC DDD parallelization note).

STEP 0 -- PAIR SMOKE (cheap, gates everything; ~minutes each, NO forest):
  for PAIR in coll_opp coll_same glissile hirth:
    JTYPE=$PAIR PAIR_SMOKE=1 LBOX=6000 NREL=500 SEED=1 OUT=smoke/$PAIR python3 evolving_forest_mfp.py
  GATE: coll_opp.length_fraction < 0.85 (annihilates, expect ~0.71 per binary_collinear.py);
        coll_same.length_fraction > 0.92 (no annihilation); glissile/hirth show n_junction>=1.
  If coll_opp does NOT annihilate, fix nf/sense BEFORE any population run (do not proceed).

STEP 1 -- FOREST-ALONE STABILITY CONTROL (one per pair x seed; bounds the unpinned-forest risk):
  for PAIR in coll_opp coll_same glissile hirth: for SEED in 1 2 3:
    JTYPE=$PAIR FOREST_ALONE=1 RHO_F=3e12 LBOX=20000 NSTEPS=2000 SEED=$SEED OUT=falone/${PAIR}_s$SEED python3 ...
  GATE: |forest_drift_frac| < 0.20 over the window. If it fails, rerun the pilot with ANCHOR=1
        (one light face-anchor per line) and re-check; if still failing, report rho_forest(t) and
        restrict the fit window to where forest is within 20% of settled.

STEP 2 -- PILOT MATRIX (smallest that can pass G1/G2):  4 pairs x 3 seeds = 12 runs, single rho_f, single tau.
  PAIRS = {coll_opp [main], coll_same [G1/G5 mechanism control], glissile, hirth}
  SEED in {1,2,3};  RHO_F=3e12 (NFOREST auto ~ derived);  TAU_MPA=45;  LBOX=20000; NGRID=64; NSTEPS=6000; REC=20.
    JTYPE=$PAIR RHO_F=3e12 TAU_MPA=45 LBOX=20000 NGRID=64 NSTEPS=6000 SEED=$SEED \
      OUT=pilot/${PAIR}_s$SEED python3 evolving_forest_mfp.py
  READOUT per run: S_store, S_remove (separate), L_mf, Lambda, R2_store, R2_remove, span_wrap_flagged,
                   mean_ambiguous_frac, gamma_global_consistency, mobile_survival trace.
  PILOT PASS: mean Lambda[coll_opp] < 0.5*mean Lambda[glissile] AND < 0.5*mean Lambda[hirth]
              AND Lambda[coll_opp] << Lambda[coll_same]  (G1+G2), with R2_store(or R2_remove)>0.8.

STEP 3 -- ONLY IF PILOT PASSES, robustness + density lever + controls (~24 runs):
  (a) rho_f lever (sqrt collapse, G2 robustness): coll_opp & glissile x RHO_F in {1e12, 1e13} x 2 seeds.
  (b) stress-consistency (G4): coll_opp & glissile x TAU_MPA in {30, 60} x 2 seeds; expect smaller Lambda
      <-> faster rho_total rise / higher needed flow stress in same run.
  (c) box + capture convergence: coll_opp & glissile x LBOX in {15000, 30000} (NGRID 32/128 to match grid
      spacing ~ maxseg) x rann in {6, 14}, 1 seed each.
  (d) G5 controls (coll_opp only, 3 seeds each, SEPARATE processes):
        FORCE=LINE_TENSION   (no pairwise elastic -> MFP shortening must weaken)
        TOPO=0               (no junction formation)
        COLL=0               (no annihilation; keep TOPO=1 so multiplication survives -> isolates the
                              annihilation channel, and is also the baseline for rho_removed cross-check)
      and coll_same is already the in-matrix sense control.
  Each G5 control raises Lambda toward the glissile/hirth value if the effect is collective-elastic.

POST-STEP ANALYSIS RECIPE (aggregating ledger.json across runs; in-script analyze() already does per-run):
  1. For each run, fit over the carrier-survival window (0.3<survival<0.95, monotone gamma):
       rho_stored(gamma_m) = c + S_store*gamma_m  ;  rho_removed_events(gamma_m) = c + S_remove*gamma_m.
  2. Baseline-subtract removal: R_remove = S_remove(opp) - S_remove(same) at matched gamma (removes
     bow/remesh/multiplication common-mode). R_store = S_store as-is (junction+residual is sense-specific).
  3. S_mf = R_store + R_remove ; L_mf = 1/(b*S_mf) ; Lambda = L_mf*sqrt(rho_f_fit).
  4. Mean Lambda per pair over seeds with CI; apply gates G1-G5.
COST: ~0.03 s/step (Trapezoid, PROGRESS); 6000 steps ~ 3 min single-thread-ish; trivially parallel.

## First sanity checks
- PAIR_SMOKE first: coll_opp must print length_fraction < 0.85 (target ~0.71 from binary_collinear.py) and coll_same > 0.92. If this fails the whole assay is invalid -- stop and fix nf/sense.
- Print list(export_data()['segs'].keys()) on the very first snapshot and confirm 'burgers' and 'planes' (or a raw matrix with >=8 cols) are present. If neither exists, the plane classifier degrades to degree-only and mean_ambiguous_frac will be large -- the run is untrustworthy.
- In the first ~200 load steps check d(gamma_m)/d(step) > 0 (mobile-only swept area is rising) AND mobile_survival starts near 1.0 then declines for coll_opp but stays ~flat for coll_same. A flat gamma_m means jamming (false infinite MFP); a survival cliff to 0 means carrier starvation (the documented free_probe failure) -- neither is a valid fit window.
- Confirm span_wrap_flagged is False over the fit window (mobile arm span < 0.45*LBOX). If True early, shrink Lseg (already 0.25*LBOX) or raise LBOX before trusting any slope.
- Check ledger CLOSURE: rho_total from the ledger vs dislocation_density(net) should agree within a few percent, and mean_ambiguous_frac < ~0.1. Large offsys/ambiguous length means the plane/Burgers tolerances are mis-binning.
- Check gamma_global_consistency = gamma_m / gamma_m_global is in a sane band (mobile-only <= global; expect 0.3-1.0). If gamma_m >> gamma_m_global the mobile filter is leaking; if << it, forest is carrying most of the strain (contamination) -- distrust L_mf.
- FOREST_ALONE control: |forest_drift_frac| < 0.20. If the unpinned forest self-coarsens/annihilates more than that without any mobile system, rho_f_fit is uninterpretable -- switch on ANCHOR=1 or restrict the window.

## Residual risks (only a real run resolves)
- export_data()['segs'] may not expose 'burgers'/'planes' on this build (no prior script reads them -- only 'nodeids' is verified used). The accessor falls back to a raw seg matrix, but if THAT is also absent the plane classifier collapses to degree-only and the collinear mobile/forest split is unrecoverable. Only a real run resolves which keys exist; first sanity check prints them.
- The step_begin/step_end hook order relative to integrate/collision/topology/remesh in THIS SimulateNetwork build is assumed (cache old positions at step_begin, swept area at step_end before surgery folds into old positions). If step_begin fires AFTER topology, the cached 'old' positions already include node surgery and the swept area is corrupted. Must verify the driver calls step_begin at the true step start; if not, move the cache into step_update_response/post_integrate. SimulateNetwork may also not define step_begin -- if so it silently never caches and gamma_m stays ~0 (the first sanity check catches this).
- get_velocities() may not exist on this ExaDisNet build; the velocity tiebreak for collinear screw/off-plane arms then never fires and those lengths default to 'mobile', possibly over-counting rho_mobile and under-counting residual. Degrades the collinear stored/mobile split, not the non-collinear pairs.
- Whether a LOW fixed tau (45 MPa) sustains FLOW without either (a) jamming on first forest contact or (b) collinear-annihilating the carriers to zero (the documented free_probe/fr_source failures) is unknown until run. The fit gates (survival 0.3-0.95, dgamma/dstep>0) reject both, but if MOST runs are rejected the sample is biased -- report the rejection fraction; if >30%, the assay is sampling-biased and tau/KM/rho_f must be retuned, or move to strain-rate control (build_pair.py idiom) which guarantees flow.
- rho_removed_events attributes ALL collision-coincident length drops to annihilation; junction zip-up of glissile/Hirth also drops length into anchored residual and could be miscounted as removed. The same-sense subtraction removes the common-mode but NOT type-specific junction zip -- so for glissile/Hirth, rho_removed_events may be slightly inflated. The stored-vs-removed separation plus FAM_J junction accounting bounds this but only a run shows the magnitude.
- The unpinned EVOLVABLE forest is the central scientific bet and the riskiest departure from all validated scripts (every prior run PINNED the forest; the one initial-density run that didn't gave KM steady state with no density lever). The FOREST_ALONE control + ANCHOR fallback + instantaneous rho_f are the mitigations, but whether an evolvable forest stays a usable, density-stable forest for long enough to fit a storage slope -- rather than collapsing to a common steady state -- can only be settled by the real run.
- NFOREST<->RHO_F mapping (NFOREST ~ RHO_F*(LBOX*b)^2) is a rough single-wrapping-line estimate; actual per-line length from insert_infinite_line depends on theta/commensurability, so the realized rho_f will differ from target. Mitigated by using the MEASURED rho_f_settled/rho_f_fit in Lambda, but the rho_f lever spacing in Step 3 may come out uneven.
- insert_infinite_line 'trial=True' kwarg is assumed (build_pair.py calls it without trial); if the signature differs the try/except defaults to inserting unconditionally, which may place a non-commensurate (open-ended) forest line. Verify the helper signature on this build.