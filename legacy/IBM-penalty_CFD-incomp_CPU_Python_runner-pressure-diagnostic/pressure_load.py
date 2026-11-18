"""Raw pressure-only surface diagnostic; no exterior-union/viscous qualification."""
import numpy as np

def pressure_load(vertices, triangles, pressure_pa):
    vertices = np.asarray(vertices, dtype=float)
    triangles = np.asarray(triangles, dtype=int)
    pressure_pa = np.asarray(pressure_pa, dtype=float)
    if vertices.ndim != 2 or vertices.shape[1] != 3 or triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError("Expected XYZ vertices and triangle connectivity")
    if triangles.min() < 0 or triangles.max() >= len(vertices):
        raise ValueError("Triangle index outside vertex array")
    if pressure_pa.shape != (len(triangles),) or not np.isfinite(vertices).all() or not np.isfinite(pressure_pa).all():
        raise ValueError("One finite original pressure value required per raw triangle")
    xyz = vertices[triangles]
    outward_area = 0.5 * np.cross(xyz[:,1]-xyz[:,0], xyz[:,2]-xyz[:,0])
    # Force on body is -p*n; +sum(p*n_x*A) is opposing +X.
    opposing_x_n = float(np.sum(pressure_pa * outward_area[:,0]))
    return {"opposing_x_N":opposing_x_n, "total_drag_qualified":False, "viscous_included":False, "exterior_union_filtered":False}
