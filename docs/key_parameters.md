# Key Tweaking Parameters

Quick reference for inversion parameters you may need to adjust.
Not exhaustive — just the ones that are easy to forget or non-obvious.

## Axis Scaling (runner.cpp setup())

| Axis | Scale | Meaning |
|------|-------|---------|
| X | 10.0 | 1 optimizer unit = 10m lateral |
| Y | 10.0 | 1 optimizer unit = 10m lateral |
| Z | 10.0 | 1 optimizer unit = 10m vertical |

Previously 50/50/10. The 5:1 XY:Z ratio made lateral movement cheaper in parameter space
to compensate for gravity's natural Z-sensitivity. Uniform 10/10/10 lets the physical
gradient drive the direction without artificial amplification.

Location: `InversionRunner::setup()` — hardcoded in the for-loop over group meshes.

## Regularization Weight (lambda)

`config.lambda` — multiplies the Laplacian smoothness penalty. Higher = stiffer surfaces.

The smoothness computes `φ = ½ Σ ||v − mean(neighbors)||²` for each free vertex.
Gradient is `λ · (v − mean(neighbors))` for own term, `−λ · (v − mean) / |N|` for neighbors.

Note: Laplacian smoothness penalizes LOCAL roughness only. Uniform expansion of the
entire mesh has zero Laplacian (all neighbors move together). To prevent global drift,
use depth bounds or reference model regularization.

## Depth Bounds (enableDepthBounds / depthBoundMargin)

```
config.enableDepthBounds = true;
config.depthBoundMargin = 2.0;  // ±20m from starting position (with scale=10)
```

Constrains every DOF to `[start − margin, start + margin]` in parameter units.
Scale-aware: margin=2.0 with axisScale=10 means ±20m physical.
This is the only mechanism that directly prevents runaway vertex migration.

## Finite-Difference Step (fdStep)

`config.fdStep = 1.0` — used only by the default `ForwardModel::computeJacobian()`.
The analytical gravity Jacobian (`GravityForward::computeGravityAnalyticJacobian`)
uses per-face FD with h=1e-4 in physical coordinates, ignoring fdStep entirely.

## L-BFGS-B History (lbfgsHistory)

`config.lbfgsHistory = 10` — number of saved gradient/position pairs for the
limited-memory Hessian approximation. Larger = better approximation but more memory.

## Armijo Constant (armijoC1)

`config.armijoC1 = 1e-4` — sufficient-decrease threshold for line search.
Smaller = stricter line search (requires more objective reduction per step).
Increase (1e-3 or 1e-2) for noisy objectives where line search rejects good steps.

## Line Search Max Iterations (lineSearchMaxIter)

`config.lineSearchMaxIter = 50` — max backtracking iterations per line search.
If the optimizer reports line search failures, increase this.

## Property Inversion (enablePropertyInversion)

When true, alternates between geometry optimization and property (density/susc/cond)
inversion. Controlled by:
- `propertyInversionInterval = 50` — geometry iters between property steps
- `propertyInversionMaxIter = 20` — max property-phase iterations
- `propertyDamping = 0.01` — Tikhonov damping for property changes

## Reference Model Regularization (enableReferenceModel / lambdaRef)

```
config.enableReferenceModel = true;
config.lambdaRef = 0.1;
```

Penalizes `½ · λ_ref · ||p − p_start||²` — deviation from starting parameters.
Unlike Laplacian smoothness, this penalizes ALL movement including uniform shifts.
Use to prevent drift while still allowing boundaries to migrate.

## Control-Point Downsampling (controlPointStride)

`config.controlPointStride = 0` — when > 0, only every N-th grid vertex gets a free DOF.
Non-control vertices are bilinearly interpolated. Reduces parameter count for speed.
0 = use all vertices.

## Mesh Boundary Mode (meshBoundaryMode)

`config.meshBoundaryMode = "free"` — controls whether hull-edge vertices are FIXED.
"free" = all vertices get vertexFreedom. "fixed" = boundary vertices are locked.

## Vertex Freedom (vertexFreedom)

`config.vertexFreedom = VertexFreedom::Z_ONLY` — default for contact surfaces.
XYZ_FREE allows full 3D movement. Z_ONLY restricts to vertical only.

This field is read by model setup (`finalizeModelSetup`) but the runner itself
never reads it — per-mesh freedoms are set during setup and then modified by
`fixExteriorFaces()`.

## Tolerance (tolerance)

`config.tolerance = 1e-6` — L-BFGS-B convergence threshold on the infinity norm
of the projected gradient. Lower = tighter convergence, more iterations.

