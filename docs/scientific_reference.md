# LithoInvert3D — Scientific Reference

**Level 3 documentation** — complete mathematical formulation of all forward models,
inversion objectives, and regularization schemes.

## 1. Joint Inversion Formulation

### 1.1 Objective Function

The inversion minimizes a weighted sum of data misfit, regularization, and constraint terms:

Φ(m) = Φ_d,grav + α_mag·Φ_d,mag + α_aem·Φ_d,aem + α_mt·Φ_d,mt + λ·Φ_r + ω·Φ_c + λ_ref·Φ_ref

where:
- Φ_d,* are per-data-type misfits (always: gravity; optional: magnetic, active EM, MT)
- Φ_r is the Laplacian surface smoothness regularization
- Φ_c is the borehole constraint penalty
- Φ_ref is the reference model deviation penalty
- α, λ, ω, λ_ref are user-set weights

### 1.2 Data Misfit (per data type)

Φ_d = ½ Σᵢ wᵢ² · (dᵢ_obs − dᵢ_pred(m))²

where wᵢ = 1/σᵢ when σᵢ > 0 (weighted), or wᵢ = 1.0 (unweighted).

### 1.3 Gradient

∇Φ = ∇Φ_d,grav + Σ α·Jᵀ·W²·r + ∇Φ_r + ∇Φ_c + ∇Φ_ref

where J is the Jacobian, W = diag(wᵢ), and r = d_obs − d_pred.

---

## 2. Forward Models

### 2.1 Gravity — Nagy Polyhedron Formula

The vertical gravity at observation point P from a homogeneous polyhedron
of density ρ is computed via surface integration over all faces:

g_z(P) = G·ρ · Σ_faces n̂_z · [ n̂·P·Ω − Σ_edges (n̂×ê)·P_edge·L ]

where:
- G = 6.67430×10⁻¹¹ m³/(kg·s²)
- n̂ = unit outward face normal, n̂_z = z-component
- Ω = solid angle subtended by the face at P
- L = ln((s₂+d₂)/(s₁+d₁)) = line integral along each edge

**Unit conversion to mGal**: g_z(mGal) = G_SI · ρ(kg/m³) · (geometric factor) × 10⁵

**Polyhedron construction**: Each litho group is bounded by:
1. Top face: upper bounding surface (CCW winding, normal up)
2. Bottom face: lower bounding surface (reversed winding, normal down)
3. Side walls: quadrilateral strips connecting boundary edges of top/bottom

### 2.2 Magnetics — Okabe/Plouff Polyhedron Formula

Total-field magnetic anomaly for a uniformly magnetized polyhedron.
Earth field direction: F̂ = (cos I·cos D, cos I·sin D, sin I).
Magnetization: M = χ·F·F̂ (induced only) or M = χ·F·F̂ + M_rem (with remanence).

Uses the same solidAngle() and lineIntegralTerm() primitives as gravity.

### 2.3 Geometry Primitives (modules/litho-core)

**Solid angle** (van Oosterom & Strackee 1983):
Ω = 2·atan2( r_a·(r_b×r_c), |r_a||r_b||r_c| + (r_a·r_b)|r_c| + (r_b·r_c)|r_a| + (r_c·r_a)|r_b| )
where r_a, r_b, r_c are vectors from observer to triangle vertices.

**Line integral term** along edge (a,b) from observer o:
L(a,b,o) = ln((s₂+d₂)/(s₁+d₁))
where s₁ = (a−o)·ê, s₂ = (b−o)·ê, d₁ = |a−o|, d₂ = |b−o|, ê = (b−a)/|b−a|.

### 2.4 EM — Status

Active EM and MT forward models currently use fallback approximations.
The full integral-equation (IE) solver is stubbed in modules/litho-em.

- **Active EM fallback**: conductive-dipole approx: response ∼ σ·δ³/R³
- **MT fallback**: 1D half-space impedance: Z = √(iωμ₀/σ)
- **Skin depth**: δ = 503/√(σ·f) [m, S/m, Hz]; time-domain: δ ≈ 503·√(t_gate/σ)

---

## 3. Regularization

### 3.1 Surface Smoothness (Laplacian Curvature)

Φ_r = ½ Σ_v ||pos(v) − mean(neighbors(v))||²

Penalizes vertices that deviate from the average position of their grid neighbors.
Flat surfaces have zero penalty; folds and creases are penalized quadratically.

### 3.2 Reference Model

Φ_ref = ½·λ_ref · ||m − m_start||²

Penalizes deviation from the starting model in parameter space.
Controlled by `[regularization]` INI section.

---

## 4. Optimization — L-BFGS-B

Limited-memory BFGS with bound constraints:

1. **Two-loop recursion** computes search direction p = −H·g without forming Hessian
2. **Line search** with backtracking Armijo + bound projection
3. **Gradient projection** for constraints: proj(x − g) − x
4. **History**: m pairs of (s_k, y_k) where s = x_{k+1}−x_k, y = g_{k+1}−g_k

Convergence criterion: ||proj(x − g) − x||_∞ < tolerance (per Byrd et al. 1995).

---

## 5. Property Inversion

Between geometry optimization phases, physical properties are optimized:

minimize ||W·(d_obs − U·ρ − U_pad·ρ_pad)||²  subject to ρ_min ≤ ρ ≤ ρ_max

where U is the unit-response matrix (response per unit property per group).
This is a bounded least-squares problem, solved via the same L-BFGS-B optimizer.

---

## 6. Durbin-Watson Diagnostic

DW = Σ(rᵢ − rᵢ₋₁)² / Σ(rᵢ)²

- DW ≈ 2.0: random residuals (good fit)
- DW ≈ 0.0: positive autocorrelation (underfitting — missing structure)
- DW ≈ 4.0: negative autocorrelation (oscillating residuals)

Two statistics computed: DW_x (sorted by x, then y) and DW_y (sorted by y, then x).

---

## 7. Terrain Correction (Bouguer Slab)

When TopographyMode::TerrainCorrected is active:
Δg_Bouguer = 2πG·ρ_bouguer · (h_station − h_datum)

where 2πG = 0.04191 mGal/(m·g/cm³).

The Bouguer slab correction is subtracted from observed gravity before misfit computation.

---

## 8. References

- Nagy, D. (1966). The gravitational attraction of a right rectangular prism. Geophysics.
- Nagy, D., Papp, G., & Benedek, J. (2000). The gravitational potential and its derivatives for the prism. J. Geodesy.
- Okabe, M. (1979). Analytical expressions for gravity anomalies due to polyhedral bodies. J. Geodetic Soc. Japan.
- Plouff, D. (1976). Gravity and magnetic fields of polygonal prisms. Geophysics.
- van Oosterom, A. & Strackee, J. (1983). The solid angle of a plane triangle. IEEE Trans. Biomed. Eng.
- Byrd, R.H., Lu, P., Nocedal, J., & Zhu, C. (1995). A limited memory algorithm for bound constrained optimization. SIAM J. Sci. Comput.
