# batch_render.py
# Run inside PyMOL:  run batch_render.py
# Or from terminal (headless): pymol -cq batch_render.py

from pymol import cmd
import os, glob, sys
from os.path import splitext, join

# -------------------------
# User-adjustable parameters
# -------------------------
RAY_WIDTH = 2400
RAY_HEIGHT = 1800
PNG_DPI = 600
OUTPUT_ROOT = "."   # working directory; change if you want outputs elsewhere
# -------------------------

def safe_name(s):
    """Make a safe PyMOL object name from a filename base (remove spaces, dots)."""
    return s.replace(" ", "_").replace(".", "_").replace("-", "_")

def render_pair(ref_obj, mob_obj, outpath):
    """Apply visuals, align, ray and save PNG to outpath (full path)."""
    try:
        # perform structural superposition (align then super for refinement)
        cmd.align(mob_obj, ref_obj)   # optional but useful as a first pass
        cmd.super(mob_obj, ref_obj)   # refinement

        # visualization
        cmd.hide("everything", "all")
        cmd.show("cartoon", ref_obj)
        cmd.show("cartoon", mob_obj)

        # transparency & colors are assumed set globally; but ensure per-object setting:
        # (the script defines 'good_orange' and 'good_grey' below)
        cmd.set("cartoon_transparency", 0.3, ref_obj)   # as you requested
        cmd.color("good_grey", ref_obj)
        cmd.color("good_orange", mob_obj)

        cmd.bg_color("white")
        cmd.set("antialias", 2)
        cmd.set("cartoon_smooth_loops", "on")
        cmd.set("cartoon_flat_sheets", "on")

        cmd.zoom("all")
        cmd.orient()

        # ray-trace then save PNG
        try:
            cmd.ray(RAY_WIDTH, RAY_HEIGHT)
            cmd.png(outpath, dpi=PNG_DPI)
            print(f"[OK] Saved (ray) {outpath}")
        except Exception as e_ray:
            # fallback: if ray fails (e.g. headless/display issues), save without ray at lower-quality
            print(f"[WARN] ray failed for {outpath}: {e_ray}. Saving fallback PNG (no ray).")
            cmd.png(outpath, dpi=PNG_DPI)
            print(f"[OK] Saved (no-ray) {outpath}")

    except Exception as e:
        print(f"[ERROR] Failed rendering pair {ref_obj} <-- {mob_obj}: {e}")

def main():
    # clean start
    cmd.reinitialize()   # clears any loaded objects

    # define colors used in the visuals (RGB 0-255)
    cmd.set_color("good_orange", [239,143,1])
    cmd.set_color("good_grey", [102,102,102])

    # find files
    ent_files = sorted(glob.glob("*.ent"))
    pdb_files = sorted(glob.glob("*.pdb"))

    if not ent_files:
        print("No .ent files found in current directory. Exiting.")
        return
    if not pdb_files:
        print("No .pdb files found in current directory. Exiting.")
        return

    print(f"Found {len(ent_files)} .ent files and {len(pdb_files)} .pdb files.")
    total_expected = len(ent_files) * len(pdb_files)
    print(f"Will produce up to {total_expected} PNG images (one per .ent x .pdb pair).")

    for i, ent_path in enumerate(ent_files, start=1):
        ent_base = splitext(os.path.basename(ent_path))[0]
        out_dir = join(OUTPUT_ROOT, ent_base)
        os.makedirs(out_dir, exist_ok=True)

        ref_obj = f"prot1_{safe_name(ent_base)}"
        print(f"\n[{i}/{len(ent_files)}] Loading reference: {ent_path} -> object '{ref_obj}'")
        try:
            cmd.load(ent_path, ref_obj)
        except Exception as e:
            print(f"[ERROR] Could not load {ent_path}: {e}")
            continue

        # iterate over pdb files
        for j, pdb_path in enumerate(pdb_files, start=1):
            pdb_base = splitext(os.path.basename(pdb_path))[0]
            mob_obj = f"prot2_{safe_name(pdb_base)}"
            out_png = join(out_dir, f"{pdb_base}.png")

            # load mobile
            try:
                cmd.load(pdb_path, mob_obj)
            except Exception as e:
                print(f"[ERROR] Could not load {pdb_path}: {e}")
                # make sure to delete any partial object
                try:
                    cmd.delete(mob_obj)
                except:
                    pass
                continue

            # render and save
            render_pair(ref_obj, mob_obj, out_png)

            # delete mobile to free memory before next iteration
            try:
                cmd.delete(mob_obj)
            except:
                pass

        # finished with this reference; delete it before next to keep namespace clean
        try:
            cmd.delete(ref_obj)
        except:
            pass

    print("\nBatch rendering finished.")

if __name__ == "__main__":
    main()
