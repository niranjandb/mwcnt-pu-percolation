"""
Monte Carlo percolation simulation for MWCNT/PU nanocomposites.
Validates 3D (bulk) vs 2D (sandwich) percolation thresholds.
CNTs modelled as randomly placed, randomly oriented capped cylinders (sticks).
Two sticks are 'connected' if their minimum surface-surface distance <= tunnelling cutoff.
Percolation = a connected cluster spans the box (bottom face to top face).
"""
import numpy as np
from scipy.spatial import cKDTree
from collections import deque

rng = np.random.default_rng(42)

# ---- CNT geometry (from TEM) ----
L = 1000.0      # nm, nanotube length
D = 20.0        # nm, outer diameter
R = D/2.0
# Tunnelling connection cutoff: two tubes conduct if gap <= delta.
# Physically ~ a few nm (tunnelling range); we use surface-surface gap cutoff.
DELTA = 3.0     # nm  (sensitivity tested later)

def random_sticks(n, box, dim):
    """Generate n sticks: centre points + unit orientation. box=(Lx,Ly,Lz)."""
    centres = rng.uniform(0, 1, size=(n, 3)) * box
    if dim == 3:
        # random 3D orientation (uniform on sphere)
        u = rng.normal(size=(n,3)); u /= np.linalg.norm(u,axis=1,keepdims=True)
    else:
        # 2D: orientation in plane (x,y), z fixed (thin slab)
        th = rng.uniform(0, np.pi, size=n)
        u = np.column_stack([np.cos(th), np.sin(th), np.zeros(n)])
        centres[:,2] = box[2]*0.5  # all in mid-plane of thin slab
    return centres, u

def seg_seg_dist(p1,d1,p2,d2,half):
    """Min distance between two line segments centred p, dir d, half-length half. Vectorized over pairs."""
    # p1,d1: (k,3) ; p2,d2:(k,3)
    r = p1 - p2
    a = np.einsum('ij,ij->i', d1, d1)*half*half*4  # not used directly; do standard clamped solve
    # Use parametric: seg1 = p1 + s*d1*half (s in [-1,1]); seg2 = p2 + t*d2*half
    A = half*d1; B = half*d2
    aa = np.einsum('ij,ij->i', A, A)
    bb = np.einsum('ij,ij->i', B, B)
    ab = np.einsum('ij,ij->i', A, B)
    ar = np.einsum('ij,ij->i', A, r)
    br = np.einsum('ij,ij->i', B, r)
    denom = aa*bb - ab*ab
    s = np.zeros_like(aa); t = np.zeros_like(aa)
    nz = denom > 1e-9
    s[nz] = (ab[nz]*br[nz] - bb[nz]*ar[nz]) / denom[nz]
    s = np.clip(s, -1, 1)
    t = (ab*s + br) / np.where(bb>1e-9, bb, 1.0)
    t = np.clip(t, -1, 1)
    s = (ab*t - ar) / np.where(aa>1e-9, aa, 1.0)
    s = np.clip(s, -1, 1)
    closest1 = p1 + (s[:,None])*A
    closest2 = p2 + (t[:,None])*B
    return np.linalg.norm(closest1-closest2, axis=1)

def build_graph(centres, u, half, cutoff_centre, gap_cutoff):
    """Connect sticks whose surface-surface gap <= gap_cutoff. Use KDTree on centres to prefilter."""
    tree = cKDTree(centres)
    pairs = tree.query_pairs(r=cutoff_centre, output_type='ndarray')
    if len(pairs)==0:
        return [[] for _ in range(len(centres))], pairs, np.array([])
    p1 = centres[pairs[:,0]]; d1 = u[pairs[:,0]]
    p2 = centres[pairs[:,1]]; d2 = u[pairs[:,1]]
    dist = seg_seg_dist(p1,d1,p2,d2,half)
    gap = dist - D  # surface-surface gap (subtract both radii = D)
    conn = gap <= gap_cutoff
    adj = [[] for _ in range(len(centres))]
    for (i,j) in pairs[conn]:
        adj[i].append(j); adj[j].append(i)
    return adj, pairs, gap

def spans(adj, centres, box, axis=2, tol=None):
    """Does a cluster touch both the low and high faces along `axis`?"""
    if tol is None: tol = L*0.5
    lo = centres[:,axis] <= (centres[:,axis].min() + tol)
    hi = centres[:,axis] >= (centres[:,axis].max() - tol)
    lo_set = set(np.where(lo)[0]); hi_set = set(np.where(hi)[0])
    seen = [False]*len(centres)
    for start in lo_set:
        if seen[start]: continue
        q=deque([start]); seen[start]=True; touch_hi=False; comp=[start]
        while q:
            n=q.popleft()
            if n in hi_set: touch_hi=True
            for m in adj[n]:
                if not seen[m]:
                    seen[m]=True; q.append(m); comp.append(m)
        if touch_hi: return True
    return False

print("Setup OK. Geometry: L=%.0f nm, D=%.0f nm, gap cutoff=%.1f nm" % (L,D,DELTA))
print("seg-seg distance + KDTree graph + spanning check defined.")

def percolation_probability(n, box, dim, trials=20):
    half = L/2.0
    cutoff_centre = L + DELTA + 1.0   # centres farther than this can't have segments within gap
    hits=0
    for _ in range(trials):
        c,u = random_sticks(n, box, dim)
        adj,_,_ = build_graph(c,u,half,cutoff_centre,DELTA)
        if spans(adj,c,box,axis=2): hits+=1
    return hits/trials

def number_to_wt(n, vol_nm3, rho_cnt=2.0, rho_pu=1.10):
    """Convert N tubes in volume (nm^3) to wt% MWCNT. Tube vol = pi r^2 L (solid cyl approx)."""
    v_tube = np.pi*(R**2)*L          # nm^3 per tube
    v_cnt = n*v_tube
    v_pu = vol_nm3 - v_cnt
    m_cnt = v_cnt*rho_cnt
    m_pu = v_pu*rho_pu
    return 100.0*m_cnt/(m_cnt+m_pu), 100.0*v_cnt/vol_nm3  # wt%, vol%

# ---------- 3D sweep ----------
print("\n=== 3D bulk percolation ===")
box3 = np.array([3000.,3000.,3000.])   # nm cube (3 x tube length)
V3 = np.prod(box3)
for n in [40,60,80,100,120,140,160]:
    p = percolation_probability(n, box3, dim=3, trials=15)
    wt,vol = number_to_wt(n, V3)
    print(f"  N={n:4d}  wt={wt:5.2f}%  vol={vol:5.3f}%  P(percolate)={p:.2f}")

print("\n=== 3D bulk percolation (corrected N range) ===")
for n in [600,800,1000,1100,1200,1400]:
    p = percolation_probability(n, box3, dim=3, trials=12)
    wt,vol = number_to_wt(n, V3)
    print(f"  N={n:4d}  wt={wt:5.2f}%  vol={vol:5.3f}%  P={p:.2f}")
