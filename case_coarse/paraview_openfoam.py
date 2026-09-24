#!/shared/ParaView/bin/pvpython
"""
OpenFOAM -> CSV processor (single script)

Implements:
1) Spatial averaging in any combination of x/y/z (0, 1, 2, or 3 directions) using TRUE reductions:
   - 0 dirs averaged: no spatial reduction (exports full dataset via SaveData)
   - 3 dirs averaged: volume-weighted mean over whole domain (single row)
   - 2 dirs averaged: area-weighted mean over the orthogonal plane, producing a 1D profile
                      along the kept axis at CELL-CENTRE coordinates (no interpolation grid)
   - 1 dir averaged: volume-weighted mean along the averaged axis, producing a 2D field
                     on the kept axes at CELL-CENTRE coordinates (no interpolation grid)

2) Time handling:
   - If TIME_AVERAGING=True: time-average ONLY over START_TIME..END_TIME, then spatial average (if enabled),
     then write ONE output CSV.
   - If TIME_AVERAGING=False: spatial average (if enabled) computed per time step, then export ALL timesteps
     in START_TIME..END_TIME, using TIMESTEP_OUTPUT_MODE: combined/separate/both.

Notes / assumptions:
- True “2 dirs averaged” (3D->1D) is computed by slicing at each unique kept-coordinate from cell centres,
  using CrinkleSlice + CellSize(Area) for area-weighting (no arbitrary sampling like Y_PROFILE_SAMPLES).
- True “1 dir averaged” (3D->2D) is computed by grouping cells by the two kept cell-centre coordinates and
  using Volume-weighted means (robust without defining a slab thickness).
"""

from paraview.simple import *
from paraview import servermanager
import os, sys, csv

# =============================================================================
# CONFIGURATION
# =============================================================================

CASE_PATH = "./"
FOAM_FILE = "foam.foam"

VARIABLES = ['U','UMean','UPrime2Mean','U_0','nut','p','pMean','pPrime2Mean','phi','phi_0']  # must exist in your case (cell arrays preferred)

# Spatial averaging configuration
# True = average in that direction (remove that dimension)
# Examples:
#  {'x': True,  'y': False, 'z': True } -> 1D profile along Y (XZ area-weighted mean at each Y cell-centre)
#  {'x': True,  'y': True,  'z': True } -> single volume-weighted mean
#  {'x': False, 'y': True,  'z': False} -> 2D field on XZ (volume-weighted mean along Y)
SPATIAL_AVERAGING = {'x': True, 'y': False, 'z': True}

# Time averaging
TIME_AVERAGING = True
START_TIME = 100.0
END_TIME = 400.0

# Output
OUTPUT_DIRNAME = "z.openfoam_output"
OUTPUT_FILENAME = "openfoam_averaged.csv"
TIMESTEP_OUTPUT_MODE = "combined"  # combined | separate | both

FIELD_PRECISION = 8

# Cell-centre grouping tolerance (mesh units). Increase if your coordinates are slightly noisy.
GROUP_TOL = 1e-9

# =============================================================================
# Utilities
# =============================================================================
def as_single_dataset(src):
    """Ensure we have a non-multiblock dataset for downstream operations."""
    mb = MergeBlocks(Input=src)
    mb.UpdatePipeline()
    return mb

def ensure_foam_file(case_path, foam_file):
    foam_path = os.path.join(case_path, foam_file)
    if not os.path.exists(foam_path):
        os.makedirs(os.path.dirname(foam_path) or ".", exist_ok=True)
        open(foam_path, "w").close()
    return foam_path

def round_to_tol(v, tol):
    if tol is None or tol <= 0:
        return v
    return round(v / tol) * tol

def write_rows_to_csv(rows, path, precision=8):
    if not rows:
        raise RuntimeError("No rows to write")

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    keys, seen = [], set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                seen.add(k); keys.append(k)

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            out = {}
            for k in keys:
                v = r.get(k, "")
                if isinstance(v, float):
                    out[k] = f"{v:.{precision}f}"
                else:
                    out[k] = v
            w.writerow(out)

def combine_csvs_with_time(separate_paths, combined_path):
    """
    separate_paths: list of (time, filepath) with identical schema.
    Writes combined CSV with a leading Time column.
    """
    os.makedirs(os.path.dirname(combined_path) or ".", exist_ok=True)

    header_written = False
    out_fields = None
    out_w = None

    with open(combined_path, "w", newline="") as out_f:
        for t, p in separate_paths:
            with open(p, "r", newline="") as in_f:
                r = csv.DictReader(in_f)
                in_fields = r.fieldnames or []
                fields = ["Time"] + in_fields

                if not header_written:
                    out_fields = fields
                    out_w = csv.DictWriter(out_f, fieldnames=out_fields)
                    out_w.writeheader()
                    header_written = True
                else:
                    if fields != out_fields:
                        raise RuntimeError(f"CSV schema changed between timesteps.\nFirst: {out_fields}\nNow:   {fields}")

                for row in r:
                    row_out = {"Time": f"{t:.6f}"}
                    row_out.update(row)
                    out_w.writerow(row_out)

def get_filtered_times(reader, start_time=None, end_time=None):
    all_times = list(reader.TimestepValues) if reader.TimestepValues else []
    if not all_times:
        return []

    tmin = start_time if start_time is not None else min(all_times)
    tmax = end_time if end_time is not None else max(all_times)
    return [t for t in all_times if tmin <= t <= tmax]

def base_name_without_ext(path):
    if "." in path:
        return path.rsplit(".", 1)[0]
    return path

def ext_or_csv(path):
    if "." in path:
        return path.rsplit(".", 1)[1]
    return "csv"

# =============================================================================
# Fetch / VTK data helpers
# =============================================================================

def fetch_dataset(source, time_value=None):
    if time_value is None:
        source.UpdatePipeline()
    else:
        source.UpdatePipeline(time_value)
    return servermanager.Fetch(source)

def get_cell_centres_points(source):
    cc = CellCenters(Input=source)
    cc.UpdatePipeline()
    vtkobj = servermanager.Fetch(cc)
    Delete(cc)

    if vtkobj is None or not hasattr(vtkobj, "GetPoints") or vtkobj.GetNumberOfPoints() == 0:
        raise RuntimeError("CellCenters produced no points")
    return vtkobj.GetPoints()

def get_cell_volumes_array(source):
    cs = CellSize(Input=source)
    cs.ComputeVolume = 1
    cs.ComputeArea = 0
    cs.ComputeLength = 0
    cs.ComputeVertexCount = 0
    cs.UpdatePipeline()
    vtkobj = servermanager.Fetch(cs)
    Delete(cs)
    if vtkobj is None or vtkobj.GetNumberOfCells() == 0:
        raise RuntimeError("CellSize produced no cells")
    vol = vtkobj.GetCellData().GetArray("Volume")
    if vol is None:
        raise RuntimeError("CellSize did not produce 'Volume' array")
    return vol

# =============================================================================
# Spatial averaging implementations (true reductions)
# =============================================================================

def volume_weighted_global_mean(source, variables):
    """
    3D -> 0D: volume-weighted mean over entire domain.
    Returns: [row]
    """
    data_vtk = fetch_dataset(source)
    cd = data_vtk.GetCellData()
    vol_arr = get_cell_volumes_array(source)

    want = set(variables)
    sumV = 0.0
    acc = {}  # name -> [sumcomp...]

    n_cells = data_vtk.GetNumberOfCells()
    for cid in range(n_cells):
        V = vol_arr.GetTuple1(cid)
        if V <= 0:
            continue
        sumV += V

        for ai in range(cd.GetNumberOfArrays()):
            arr = cd.GetArray(ai)
            name = arr.GetName()
            if not name:
                continue
            base = name.replace("_average", "")
            if base not in want:
                continue

            ncomp = arr.GetNumberOfComponents()
            tup = arr.GetTuple(cid)

            if name not in acc:
                acc[name] = [0.0] * ncomp
            for c in range(ncomp):
                acc[name][c] += V * tup[c]

    row = {"Volume": sumV}
    for name, sums in acc.items():
        if len(sums) == 1:
            row[f"{name}_mean"] = sums[0] / sumV if sumV else float("nan")
        else:
            for c, s in enumerate(sums):
                row[f"{name}_mean_{c}"] = s / sumV if sumV else float("nan")

    return [row]

def slice_area_weighted_mean_at_coord(source, normal_axis, coord_value, variables):
    axis_to_normal = {'x':[1,0,0], 'y':[0,1,0], 'z':[0,0,1]}
    idx = {'x':0,'y':1,'z':2}[normal_axis]
    normal = axis_to_normal[normal_axis]

    source.UpdatePipeline()
    xmin, xmax, ymin, ymax, zmin, zmax = source.GetDataInformation().GetBounds()
    origin = [(xmin+xmax)/2.0, (ymin+ymax)/2.0, (zmin+zmax)/2.0]
    origin[idx] = coord_value

    sl = Slice(Input=source)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = origin
    sl.SliceType.Normal = normal
    sl.UpdatePipeline()

    integ = IntegrateVariables(Input=sl)
    integ.UpdatePipeline()

    vtkobj = servermanager.Fetch(integ)
    Delete(integ); Delete(sl)

    if vtkobj is None or vtkobj.GetNumberOfCells() == 0:
        return None

    cd = vtkobj.GetCellData()
    area_arr = cd.GetArray("Area")
    if area_arr is None:
        # debug to see what arrays exist
        names = [cd.GetArray(i).GetName() for i in range(cd.GetNumberOfArrays())]
        # print(f"[DEBUG] IntegrateVariables arrays at {normal_axis}={coord_value}: {names}")
        return None

    A = area_arr.GetTuple1(0)
    if A == 0.0:
        return None

    want = set(variables)
    out = {"Area": A}

    for ai in range(cd.GetNumberOfArrays()):
        arr = cd.GetArray(ai)
        name = arr.GetName()
        if not name or name == "Area":
            continue

        base = name.replace("_average", "")
        if base not in want:
            continue

        ncomp = arr.GetNumberOfComponents()
        tup = arr.GetTuple(0)

        if ncomp == 1:
            out[f"{name}_mean"] = tup[0] / A
        else:
            for c in range(ncomp):
                out[f"{name}_mean_{c}"] = tup[c] / A

    return out

def profile_area_weighted_mean(source, kept_axis, variables, tol=1e-9):
    """
    3D -> 1D: two directions averaged, one kept.
    For each unique cell-centre coordinate along kept_axis, compute an area-weighted mean
    on the orthogonal plane at that coordinate.
    Returns rows with kept coordinate (X or Y or Z) + means.
    """
    pts = get_cell_centres_points(source)
    n = pts.GetNumberOfPoints()
    ax_i = {'x':0,'y':1,'z':2}[kept_axis]

    coords = set()
    for i in range(n):
        coords.add(round_to_tol(pts.GetPoint(i)[ax_i], tol))
    coords = sorted(coords)

    # print(f"[DEBUG] kept_axis={kept_axis} unique coords={len(coords)} "
    #   f"min={coords[0] if coords else None} max={coords[-1] if coords else None}")
    
    rows = []
    for c in coords:
        mean = slice_area_weighted_mean_at_coord(source, normal_axis=kept_axis, coord_value=c, variables=variables)
        if mean is None:
            continue
        row = {kept_axis.upper(): float(c)}
        row.update(mean)
        rows.append(row)

    rows.sort(key=lambda r: r[kept_axis.upper()])
    return rows

def field2d_volume_weighted_mean(source, avg_axis, variables, tol=1e-9):
    """
    3D -> 2D: one direction averaged, two kept.
    Groups cells by the two kept cell-centre coordinates and computes volume-weighted means.
    Returns rows with kept coords and averaged variables (no interpolation grid).
    """
    data_vtk = fetch_dataset(source)
    cd = data_vtk.GetCellData()
    vol_arr = get_cell_volumes_array(source)
    pts = get_cell_centres_points(source)

    want = set(variables)
    ax = {'x':0,'y':1,'z':2}
    kept_axes = [a for a in ['x','y','z'] if a != avg_axis]
    i1, i2 = ax[kept_axes[0]], ax[kept_axes[1]]

    # accumulators keyed by (c1,c2)
    sumV = {}
    acc = {}  # (c1,c2) -> {name -> [sumcomp...]}

    n_cells = data_vtk.GetNumberOfCells()
    for cid in range(n_cells):
        V = vol_arr.GetTuple1(cid)
        if V <= 0:
            continue

        p = pts.GetPoint(cid)
        c1 = round_to_tol(p[i1], tol)
        c2 = round_to_tol(p[i2], tol)
        key = (c1, c2)

        sumV[key] = sumV.get(key, 0.0) + V
        if key not in acc:
            acc[key] = {}

        for ai in range(cd.GetNumberOfArrays()):
            arr = cd.GetArray(ai)
            name = arr.GetName()
            if not name:
                continue
            base = name.replace("_average", "")
            if base not in want:
                continue

            ncomp = arr.GetNumberOfComponents()
            tup = arr.GetTuple(cid)

            if name not in acc[key]:
                acc[key][name] = [0.0] * ncomp
            for c in range(ncomp):
                acc[key][name][c] += V * tup[c]

    rows = []
    for (c1, c2), V in sumV.items():
        row = {kept_axes[0].upper(): float(c1), kept_axes[1].upper(): float(c2), "Volume": V}
        for name, sums in acc[(c1,c2)].items():
            if len(sums) == 1:
                row[f"{name}_mean"] = sums[0] / V
            else:
                for k, s in enumerate(sums):
                    row[f"{name}_mean_{k}"] = s / V
        rows.append(row)

    rows.sort(key=lambda r: (r[kept_axes[0].upper()], r[kept_axes[1].upper()]))
    return rows

def spatial_average_any(source, spatial_avg_cfg, variables, tol=1e-9):
    """
    Returns ('proxy', proxy) or ('rows', rows).

    - 0 averaged dirs: proxy = source (no reduction)
    - 1 averaged dir: rows = 2D field on kept plane (volume-weighted)
    - 2 averaged dirs: rows = 1D profile on kept axis (area-weighted)
    - 3 averaged dirs: rows = single row (volume-weighted)
    """
    avg_axes = [a for a, do in spatial_avg_cfg.items() if do]
    keep_axes = [a for a in ['x','y','z'] if a not in avg_axes]

    if len(avg_axes) == 0:
        return ("proxy", source)

    if len(avg_axes) == 3:
        return ("rows", volume_weighted_global_mean(source, variables))

    if len(avg_axes) == 2:
        kept_axis = keep_axes[0]
        return ("rows", profile_area_weighted_mean(source, kept_axis, variables, tol))

    if len(avg_axes) == 1:
        avg_axis = avg_axes[0]
        return ("rows", field2d_volume_weighted_mean(source, avg_axis, variables, tol))

    raise RuntimeError("Unexpected averaging configuration")

# =============================================================================
# Time averaging over selected range (no TemporalStatistics)
# =============================================================================

def time_average_rows_over_times(per_time_rows):
    """
    per_time_rows: list of rows (dicts) per time step for the SAME spatial keys.
    This averages numeric fields across time for each spatial key.
    Spatial key is inferred as all non-numeric columns? Too fragile.
    Here we use all coordinate-like columns among X,Y,Z (and both for 2D).
    """
    if not per_time_rows:
        return []

    # Determine which coordinate columns exist
    coord_cols = [c for c in ["X","Y","Z"] if c in per_time_rows[0]]
    if not coord_cols:
        # 0D case: no coords
        coord_cols = []

    sums = {}
    counts = {}

    for r in per_time_rows:
        key = tuple(r.get(c) for c in coord_cols) if coord_cols else ("__global__",)
        sums.setdefault(key, {})
        counts[key] = counts.get(key, 0) + 1

        for k, v in r.items():
            if k in coord_cols:
                continue
            if isinstance(v, float):
                sums[key][k] = sums[key].get(k, 0.0) + v

    out = []
    for key in sums:
        row = {}
        if coord_cols:
            for i, c in enumerate(coord_cols):
                row[c] = key[i]
        for k, s in sums[key].items():
            row[k] = s / float(counts[key])
        out.append(row)

    # stable sort by coords
    if coord_cols:
        out.sort(key=lambda r: tuple(r[c] for c in coord_cols))
    return out

# =============================================================================
# Main
# =============================================================================

def main():
    print("="*80)
    print("ParaView OpenFOAM Processor (true spatial averaging any direction)")
    print("="*80)

    foam_path = ensure_foam_file(CASE_PATH, FOAM_FILE)
    out_dir = os.path.join(CASE_PATH, OUTPUT_DIRNAME)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, OUTPUT_FILENAME)

    reader = OpenFOAMReader(FileName=foam_path)
    reader.MeshRegions = ["internalMesh"]
    reader.CellArrays = VARIABLES

    # Try to force reading of 0/ (initial conditions)
    if hasattr(reader, "SkipZeroTime"):
        reader.SkipZeroTime = 0
    if hasattr(reader, "ReadZeroTime"):
        reader.ReadZeroTime = 1

    reader.UpdatePipeline()

    time_steps = get_filtered_times(reader, START_TIME, END_TIME)
    print(f"Selected timesteps: {time_steps}")
    
    if time_steps:
        print(f"Time range selected: {min(time_steps):.6f} .. {max(time_steps):.6f}")

    print(f"SPATIAL_AVERAGING: {SPATIAL_AVERAGING}")
    print(f"TIME_AVERAGING: {TIME_AVERAGING}")
    print(f"TIMESTEP_OUTPUT_MODE: {TIMESTEP_OUTPUT_MODE}")
    print(f"Output: {out_path}")

    current = MergeBlocks(Input=reader)
    current.UpdatePipeline()

    tk = GetTimeKeeper()

    # -------------------------------------------------------------------------
    # TIME AVERAGING ENABLED: average only over START..END, then spatial average once
    # -------------------------------------------------------------------------
    if TIME_AVERAGING:
        if not time_steps:
            raise RuntimeError("No timesteps in selected START_TIME..END_TIME range")

        # Compute spatially reduced rows for each time step, then time-average those rows
        all_rows = []
        for t in time_steps:
            tk.Time = t
            reader.UpdatePipeline(t)

            kind, obj = spatial_average_any(current, SPATIAL_AVERAGING, VARIABLES, GROUP_TOL)
            if kind == "proxy":
                # No spatial reduction requested: fall back to SaveData (not row-based)
                # For time averaging in this mode, user would need a different representation.
                raise RuntimeError("TIME_AVERAGING=True with no spatial reduction is not supported in this script.")
            else:
                all_rows.append(obj)

        # Flatten and time-average by matching coordinate keys
        # (We assume each time step produces the same set of coordinate rows.)
        # Average row-by-row using coordinate columns as key.
        flat = []
        for rows in all_rows:
            flat.extend(rows)

        # Better: average per key across time
        # Build per-time list then combine:
        # We'll merge time by key explicitly:
        coord_cols = [c for c in ["X","Y","Z"] if c in all_rows[0][0]]
        if not coord_cols:
            coord_cols = []

        sums = {}
        counts = {}

        for rows in all_rows:
            for r in rows:
                key = tuple(r.get(c) for c in coord_cols) if coord_cols else ("__global__",)
                sums.setdefault(key, {})
                counts[key] = counts.get(key, 0) + 1
                # keep coords
                for c in coord_cols:
                    sums[key][c] = r.get(c)
                for k, v in r.items():
                    if k in coord_cols:
                        continue
                    if isinstance(v, float):
                        sums[key][k] = sums[key].get(k, 0.0) + v

        out_rows = []
        for key, d in sums.items():
            row = {}
            for c in coord_cols:
                row[c] = d.get(c)
            for k, v in d.items():
                if k in coord_cols:
                    continue
                if isinstance(v, float):
                    row[k] = v / float(counts[key])
            out_rows.append(row)

        if coord_cols:
            out_rows.sort(key=lambda r: tuple(r[c] for c in coord_cols))

        write_rows_to_csv(out_rows, out_path, FIELD_PRECISION)
        print("✓ Wrote time-averaged (range-limited) result.")
        return

    # -------------------------------------------------------------------------
    # TIME AVERAGING DISABLED: spatial average per time, export all selected times
    # -------------------------------------------------------------------------
    if not time_steps:
        # single/static case: just current time
        kind, obj = spatial_average_any(current, SPATIAL_AVERAGING, VARIABLES, GROUP_TOL)
        if kind == "proxy":
            SaveData(out_path, proxy=obj, Precision=FIELD_PRECISION)
        else:
            write_rows_to_csv(obj, out_path, FIELD_PRECISION)
        print("✓ Wrote single-time output.")
        return

    base = base_name_without_ext(out_path)
    ext = ext_or_csv(out_path)

    separate_paths = []

    for i, t in enumerate(time_steps):
        tk.Time = t
        current.UpdatePipeline(t)

        kind, obj = spatial_average_any(current, SPATIAL_AVERAGING, VARIABLES, GROUP_TOL)

        if kind == "proxy":
            # no spatial reduction: SaveData full dataset per time
            p = f"{base}_t{i:04d}_time{t:.6f}.{ext}"
            SaveData(p, proxy=obj, Precision=FIELD_PRECISION, WriteTimeSteps=0)
        else:
            p = f"{base}_t{i:04d}_time{t:.6f}.{ext}"
            write_rows_to_csv(obj, p, FIELD_PRECISION)

        separate_paths.append((t, p))

        if (i + 1) % 10 == 0 or i == len(time_steps) - 1:
            print(f"  Saved {i+1}/{len(time_steps)}")

    if TIMESTEP_OUTPUT_MODE in ("combined", "both"):
        combine_csvs_with_time(separate_paths, out_path)
        print(f"✓ Combined CSV saved: {out_path}")

    if TIMESTEP_OUTPUT_MODE in ("separate", "both"):
        print(f"✓ Separate CSVs saved: {base}_t*_time*.{ext}")

    if TIMESTEP_OUTPUT_MODE not in ("combined", "separate", "both"):
        raise ValueError("TIMESTEP_OUTPUT_MODE must be: combined, separate, or both")

    print("✓ Done.")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
