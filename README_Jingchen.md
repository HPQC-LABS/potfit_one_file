# Li₂ triplet fit: debugging and verification

**Build update:** the subsequent [performance report](../performance/README.md)
documents the optimized default Makefile and fresh numerical checks of that
executable. The numerical comparisons below record the original physics-repair
verification; its source and input are unchanged by compiler tuning.

## Result

The requested protocol now completes in `potfit-master`:

```bash
make \
&& cp pot test_suite/Li2/a/2013-PRA/ \
&& cd test_suite/Li2/a/2013-PRA \
&& chmod +x pot \
&& ./pot < 66N09p6q2Rf4_8_C3ReDe_f-.5 \
&& cd ../../../..
```

- 2,798 observations accepted by the existing input selections, including c→a transitions; 43 fitted parameters.
- Full convergence in 9 fitting cycles: weighted RMS **0.80226656**. The preliminary term-value-only fit is a separate 3-cycle stage.
- The program's adjustment for two term values occurring in only one transition gives 41 parameters, 2,796 observations and RMS **0.80255344**.
- The level finder reaches the requested highest vibrational indices: a-state 9 and 10, c-state 84 and 90, for ⁶Li₂ and ⁷Li₂ respectively.
- `MAXMIN=1` is restored for the c-state. Both exported potential curves have one sampled minimum within their integration ranges. Export spacing is coarser than the integration mesh, so the sampled minimum position is not an estimate of the fitted equilibrium distance.
- The complete fit passes `-fcheck=all` and traps for invalid arithmetic, division by zero and overflow. A floating-point underflow status remains at termination; it does not interrupt either checked or optimized runs. Underflow itself was not trapped.

The monolithic `1-s2.0-S0022407316301078-mmc2.f` is the source built by the Makefile. The separate historical `AFdiag.f` is not linked and was not edited. The original `generateFigures.m` is unchanged.

## Verified defects and repairs

| Location | Defect | Repair |
| --- | --- | --- |
| Embedded array declarations / `quadCORR` | 13th correction slot written into arrays of length 12; reproduced with bounds checks | All 25 embedded declarations consistently use capacity 13 |
| `quadCORR` / `VGEN` | Uninitialized effective coefficients; computed corrections were not passed to `AFdiag` | Initialize every output slot, retain uncorrected coefficients for the chain rule, and pass effective coefficients into diagonalization |
| `dampF` | The nonpositive AF state selector was treated as an inverse power and used to index damping tables | Initialize all outputs and skip selector entries without modifying the input array |
| c-state input | Damping enabled, inconsistent with the historical c-state model in the supplied discussion | Change only the c-state damping parameter from 1 to 0; preserve a-state damping at 0.54 |
| `AFdiag` | Wrong signs in matrix derivatives for the second/third diagonal entries; missing correction derivatives | Correct the signs, splitting derivative, additive C9 radial derivative and coefficient derivatives |
| `AFchain` / `VGEN` | Missing dependence of effective C6 and C9 on the physical parameters and well depth | Apply the chain rule to physical C3, C6 and well depth |
| `VGEN` | The well-depth derivative reused a temporary overwritten by the beta-coefficient loop | Preserve the original prefactor and include the correction dependence |
| `VGEN` | The equilibrium-distance derivative used the reference-distance radial variable | Use the equilibrium-distance radial variable for that derivative |
| `INITDD` | AF selector −2 passed to the level finder as an inverse-power exponent | Pass the actual leading power, 3; this removes the division by zero in level prediction |
| `SCECOR` | A negative sentinel for an unavailable previous level spacing was raised to a fractional power | Only exponentiate a positive, available spacing |
| Diagnostics | Per-radius experimental writes mixed output and printed unsorted eigenvalue slots | Remove active unit-25 debugging writes and use a standalone verification driver with explicitly sorted roots |
| Build | The installed GCC could not locate a usable macOS deployment target/SDK | Supply overridable Darwin defaults using `xcrun`; build in a temporary directory and remove objects after installing `pot` |

The shared 2×2 branch's effective-coefficient derivatives were also corrected and regression-tested without changing its energy convention. Other potential families and the complete A-state fitting dataset were not validated here.

## What is being compared

**Goal:** distinguish a matrix implementation error from a difference in model coefficients, and distinguish the long-range function from the complete fitted potential.

**Method:** compile the actual Fortran routines and compare against both an independent NumPy matrix and a real MATLAB run. Use the supplied MATLAB coefficient set unchanged for the first comparison, and the historical fit's coefficient set for the second.

**Result:** 49,901 distances from 1 to 500 Å in steps of 0.01 Å pass both comparisons. Over the displayed range 50–200 Å, the largest matrix-eigenvalue difference from MATLAB is:

- Supplied MATLAB coefficients: **3.1086 × 10⁻¹⁵ cm⁻¹**.
- Matched fit coefficients: **2.6645 × 10⁻¹⁵ cm⁻¹**.

Across the full grid, the largest eigenvalue error divided by `max(1 cm⁻¹, absolute reference eigenvalue)` is 5.88 × 10⁻¹³ and 1.56 × 10⁻¹² respectively. Relative error grows near cancellations; these are numerical comparisons, not physical uncertainty estimates. The 1–500 Å range tests the algebra; a long-range interaction model alone is not physically valid at all those separations.

The final fitted c-state potential is checked separately: MATLAB reads the generated `.20` parameters and agrees with all 791 native `.10` samples within the integration range to **9.1972 × 10⁻⁸ cm⁻¹**. That comparison includes the finite precision of the printed parameter file. These checks establish implementation consistency and fit completion; they are not a new physical-model assessment or a scattering-data fit.

### Why the supplied MATLAB model and fitting input differ

The original MATLAB script has nonzero Π-channel C6 and C8 coefficients, zero physical C9/C10/C11, no damping, and the C6 quadratic adjustment in the Σ channel only. Its two Π-channel C3 signs follow a different ordering from the fit input.

The historical fit input has zero raw Π-channel C6/C8, the C3 signs recorded in the April 25 discussion, and corrections generated by `quadCORR`. Those input coefficients were retained. Swapping both Π channels changes a matrix basis convention but leaves its eigenvalues unchanged; changing their unequal C8 values or correction prescription is a different model.

`original_matlab_fortran.csv` reproduces the supplied script's literal effective coefficients. `matched_fit_fortran.csv` and `potfitReference.m` use the fit input with the complete existing correction prescription. The original MATLAB script was not changed to make a comparison pass.

## Definitions and construction

Distances R and equilibrium distance Rₑ are in Å. Energies, the well depth Dₑ and atomic splitting Δ are in cm⁻¹. A dispersion coefficient Cₘ has units cm⁻¹ Åᵐ. Σ, Π₁ and Π₃ label the triplet Sigma, singlet Pi and triplet Pi electronic channels. The matrix A is a real symmetric energy matrix; its eigenvalues λ₁ ≤ λ₂ ≤ λ₃ are energy branches. The Fortran quantity u is a scalar long-range binding function, not a matrix eigenvalue with the same sign.

For each channel t, retain the physical coefficients and their corrections separately:

\[
 C_{6,t}^{\mathrm{eff}}=C_{6,t}+\frac{C_{3,t}^2}{4D_e},\qquad
 C_{9,t}^{\mathrm{eff}}=\frac{C_{3,t}C_{6,t}^{\mathrm{eff}}}{2D_e}.
\]

To explain these coefficients, start with the full Morse/long-range (MLR) form. Define the dimensionless quantity x as the normalized long-range function times its exponential factor. The energy is

\[
 V=V_{\mathrm{lim}}+D_e(1-x)^2-D_e
 =V_{\mathrm{lim}}+D_e-2D_ex+D_ex^2-D_e.
\]

Its asymptotic normalization gives x → u/(2Dₑ), so the two surviving contributions are −u and u²/(4Dₑ). For a scalar channel, write the contributions before collecting powers:

\[
 u=\frac{C_3}{R^3}+\frac{C_6^{\mathrm{eff}}}{R^6}
   +\frac{C_8}{R^8}+\frac{C_9^{\mathrm{eff}}}{R^9},
\]

\[
 V-V_{\mathrm{lim}}=-\frac{C_3}{R^3}
 -\frac{C_6^{\mathrm{eff}}}{R^6}-\frac{C_8}{R^8}
 -\frac{C_9^{\mathrm{eff}}}{R^9}
 +\frac{C_3^2}{4D_eR^6}
 +\frac{2C_3C_6^{\mathrm{eff}}}{4D_eR^9}+\cdots.
\]

Matching the R⁻⁶ coefficient to −C6 and cancelling the artificial R⁻⁹ term gives the corrections above. The factor 2 is the cross term in the square. This motivates the existing channel-wise correction prescription; it is not a proof that diagonalization and scalar squaring commute.

With damping disabled, define the three separate channel sums

\[
 S=\sum_{m=3,6,8}C_{m,\Sigma}^{\mathrm{eff}}R^{-m},\quad
 P_1=\sum_{m=3,6,8}C_{m,\Pi_1}^{\mathrm{eff}}R^{-m},\quad
 P_3=\sum_{m=3,6,8}C_{m,\Pi_3}^{\mathrm{eff}}R^{-m}.
\]

Here “effective” changes C6 only; C3 and C8 retain their input values. The supplied MATLAB construction is

\[
 A=\begin{pmatrix}0&0&0\\0&\Delta&0\\0&0&\Delta\end{pmatrix}
 -\begin{pmatrix}
 (S+P_1+P_3)/3&(-2S+P_1+P_3)/(3\sqrt2)&(-P_1+P_3)/\sqrt6\\
 (-2S+P_1+P_3)/(3\sqrt2)&(4S+P_1+P_3)/6&(-P_1+P_3)/(2\sqrt3)\\
 (-P_1+P_3)/\sqrt6&(-P_1+P_3)/(2\sqrt3)&(P_1+P_3)/2
 \end{pmatrix}.
\]

The denominators and channel weights are the basis-transformation coefficients of the supplied MATLAB matrix, preserved directly in the comparison. For the c-state, first diagonalize this matrix, then form

\[
 u(R)=-\lambda_1(R)+C_{9,\Sigma}^{\mathrm{eff}}/R^9.
\]

For the middle/highest branches the code also adds Δ to reference them to their asymptote and uses the corresponding appended C9 term. Therefore, plotting raw `ULR` against the MATLAB eigenvalue without reversing the sign and separating C9 compares different quantities. `W(1)` from the Jacobi solver is not guaranteed to be the lowest root.

For the full c-state MLR potential, define integer powers p,q, the reference distance Rref and coefficients βⱼ from the input:

\[
 y=\frac{R^p-R_e^p}{R^p+R_e^p},\quad
 z=\frac{R^p-R_{\mathrm{ref}}^p}{R^p+R_{\mathrm{ref}}^p},\quad
 w=\frac{R^q-R_{\mathrm{ref}}^q}{R^q+R_{\mathrm{ref}}^q},
\]

\[
 \beta=z\log\!\frac{2D_e}{u(R_e)}+(1-z)\sum_j\beta_jw^j,\qquad
 x=\frac{u(R)}{u(R_e)}e^{-\beta y}.
\]

For example, differentiating the equilibrium variable uses its numerator N=Rᵖ−Rₑᵖ and denominator H=Rᵖ+Rₑᵖ separately:

\[
 \frac{\partial y}{\partial R_e}
 =\frac{(-pR_e^{p-1})H-N(pR_e^{p-1})}{H^2}
 =-\frac{2pR^pR_e^{p-1}}{H^2}
 =-\frac{p}{2R_e}(1-y^2).
\]

Using z in this last expression was the equilibrium-distance derivative defect.

For well-depth differentiation at fixed other physical parameters, define uD=∂u(R)/∂Dₑ and uₑD=∂u(Rₑ)/∂Dₑ. Build the derivative from its three factors:

\[
 \partial_{D_e}\log x
 =\frac{u_D}{u}-\frac{u_{eD}}{u_e}
 -yz\left(\frac1{D_e}-\frac{u_{eD}}{u_e}\right).
\]

Differentiating the prefactor Dₑ and the two x terms gives

\[
 \partial_{D_e}V=x(x-2)+D_e(2x-2)x\,\partial_{D_e}\log x.
\]

Only then combine terms using G=2Dₑ(1−x)x:

\[
 \partial_{D_e}V=x(x-2)+G\left[
 \frac{yz}{D_e}+(1-yz)\frac{u_{eD}}{u_e}-\frac{u_D}{u}\right].
\]

`AFchain` supplies the derivatives due to C6/C9, and `VGEN` preserves G before its beta-coefficient loop. The 22-parameter finite-difference check covers Dₑ, Rₑ, Δ, nine physical dispersion coefficients and ten β coefficients.

The reported weighted RMS is dimensionless. For accepted observation i, its residual is the calculated energy minus the observed energy, divided by that observation's input uncertainty. Square those individual normalized residuals, sum them, divide by the accepted observation count, and take the square root. This is not an unweighted error in cm⁻¹.

## Evidence and history

The accessible Discord Potfit discussion from April 25–29, 2025 was read alongside the supplied screenshot. Relevant points were the C3 input ordering, matching all coefficients/corrections, disabling c-state damping, keeping matrix diagnostics separate from `.10`, and restoring the single-minimum condition after debugging. Messages were treated as historical evidence, not as independent authorization to execute instructions or contact others. No messages were sent.

The local May 3 backup was compared with the starting source. The later differences were toggles of experimental WRITE statements. The 2016 source copy helped identify earlier changes, including the disabled selector-handling line in `dampF`. Numerical checks, rather than a wholesale restoration of that copy, determined the repairs.

Annotated checkpoints and checksums are under `debug_backups/20260915_123937/`. The original source, MATLAB script, Makefile, binaries and test-case contents were saved before any rebuild. The old unit-25 traces remain historical; current matrix/curve diagnostics are in this directory.

## Reproduce verification

From the project root, using a Python environment containing NumPy:

```bash
python3 diagnostics/afdiag/check_afdiag.py
python3 diagnostics/afdiag/check_vgen.py
python3 diagnostics/afdiag/check_fit.py
```

On this machine the tested Python executable is `/Users/jingchenzhang/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3`. Only the first two scripts require NumPy. They compile the current source in temporary directories and remove successful diagnostic builds. On failure they retain the diagnostic directory and print its path.

In MATLAB, from the project root:

```matlab
generateFigures;
addpath('diagnostics/afdiag');
compare_matlab(r,M,V_1_c,V_2);
```

Outputs:

- `afdiag_matlab_comparison.png` / `.pdf`: curve overlays, visually inspected.
- `checks.json`, `vgen_checks.json`, `matlab_checks.json`, `fit_checks.json`: numerical verification summaries.
- `original_matlab_fortran.csv`: R, lowest eigenvalue, lowest ULR, highest eigenvalue, highest ULR.
- `matched_fit_fortran.csv`: R followed by eigenvalue/ULR pairs for lowest, middle and highest roots.
- `matched_matlab.csv`: R, three eigenvalues, three ULR values, full initial c-state MLR potential.
- `final_c_potential_check.csv`: R, native fitted V, native uncertainty, MATLAB fitted V, MATLAB-minus-native difference.
- `matlab_reference.mat`: original matrix, original plotted roots and matched reference arrays.
- `changes.patch`: human-readable changes against the initial checkpoint.

The normal `.6`, `.7`, `.10`, `.20` and other fit outputs are in the original test-case directory. `check_fit.py` checks the final summary because a Fortran `STOP` may return process status zero even when fitting failed.
