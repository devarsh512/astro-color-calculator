# =============================================================
# Astronomy Color Calculator (Final Version)
# With SDSS Flux Fixes, FITS Support, and Comparison Tools
# =============================================================

import os, io, uuid, time, math, requests # pyright: ignore[reportMissingModuleSource]
import numpy as np # pyright: ignore[reportMissingImports]
import matplotlib.pyplot as plt # pyright: ignore[reportMissingModuleSource]

from flask import Flask, render_template, request # pyright: ignore[reportMissingImports]
from astropy.io import fits # pyright: ignore[reportMissingImports]
from astropy import units as u # pyright: ignore[reportMissingImports]
from astropy.constants import c # pyright: ignore[reportMissingImports]

app = Flask(__name__)

BASE = os.path.dirname(os.path.abspath(__file__))
FILTER_DIR = os.path.join(BASE, "filters")
STATIC_DIR = os.path.join(BASE, "static")
os.makedirs(STATIC_DIR, exist_ok=True)

# -------------------------------------------------------------
# AUTO CLEAN STATIC FOLDER
# -------------------------------------------------------------
def clean_static_folder(max_age_minutes=10):
    now = time.time()
    for f in os.listdir(STATIC_DIR):
        if f.endswith(".png"):
            fp = os.path.join(STATIC_DIR, f)
            if os.path.isfile(fp):
                if (now - os.path.getmtime(fp)) > max_age_minutes * 60:
                    try:
                        os.remove(fp)
                    except:
                        pass

# -------------------------------------------------------------
# LOAD FILTERS
# -------------------------------------------------------------
def load_filter(path):
    arr = np.loadtxt(path)
    return arr[:,0] * u.AA, arr[:,1]

u_w, u_T = load_filter(os.path.join(FILTER_DIR, "sdss_u.dat"))
g_w, g_T = load_filter(os.path.join(FILTER_DIR, "sdss_g.dat"))
r_w, r_T = load_filter(os.path.join(FILTER_DIR, "sdss_r.dat"))
i_w, i_T = load_filter(os.path.join(FILTER_DIR, "sdss_i.dat"))
z_w, z_T = load_filter(os.path.join(FILTER_DIR, "sdss_z.dat"))

FILTERS = {"u":(u_w,u_T), "g":(g_w,g_T), "r":(r_w,r_T),
           "i":(i_w,i_T), "z":(z_w,z_T)}

# -------------------------------------------------------------
# SYNTHETIC AB MAGNITUDE CALCULATION
# -------------------------------------------------------------
def synthetic_abmag(w, f_lambda, fw, ft):
    try:
        sed_w = w.to(u.AA).value
        sed_f = f_lambda.to(u.erg/(u.s*u.cm**2*u.AA)).value
        fil_w = fw.to(u.AA).value

        f_interp = np.interp(fil_w, sed_w, sed_f, left=0, right=0)
        num = np.trapz(f_interp * ft * fil_w, fil_w)
        den = np.trapz(ft * fil_w, fil_w)

        if den == 0 or num <= 0:
            return None

        c_ang = c.to(u.AA/u.s).value
        f_nu = (num/den) / c_ang
        return -2.5 * np.log10(f_nu) - 48.60

    except:
        return None

# -------------------------------------------------------------
# BLACKBODY SHAPE
# -------------------------------------------------------------
def blackbody_shape(T, n=2000):
    w = np.linspace(3000, 11000, n) * u.AA
    lam_cm = w.to(u.cm).value
    expo = (1.4388e-2) / ((w.to(u.m)).value * T)
    expo = np.minimum(expo, 700)
    B = (1/lam_cm**5) / (np.exp(expo) - 1)
    return w, B

# -------------------------------------------------------------
# READ SDSS FITS SPECTRUM
# -------------------------------------------------------------
def read_sdss_fits(fileobj):
    fileobj.seek(0)
    hd = fits.open(io.BytesIO(fileobj.read()))
    data = hd[1].data
    lam = 10**data["loglam"] * u.AA
    flux_raw = data["flux"]             # SDSS flux in (1e-17 erg/s/cm^2/Å)/sr
    hd.close()
    return lam, flux_raw

# -------------------------------------------------------------
# READ ASCII TEXT SED
# -------------------------------------------------------------
def read_text_sed(fileobj):
    fileobj.seek(0)
    raw = fileobj.read()
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8", errors="ignore")
    arr = np.loadtxt(io.StringIO(raw))
    return arr[:,0]*u.AA, arr[:,1]*(u.erg/(u.s*u.cm**2*u.AA))

# -------------------------------------------------------------
# DOWNLOAD FITS VIA URL
# -------------------------------------------------------------
def download_fits_url(url):
    r = requests.get(url, timeout=15)
    if r.status_code != 200:
        raise ValueError("Download failed")
    return io.BytesIO(r.content)

# -------------------------------------------------------------
# PLOT SED + FILTERS
# -------------------------------------------------------------
def plot_sed(w, f):
    plt.figure(figsize=(7,4))
    plt.plot(w.value, f.value, lw=1.2, label="SED")

    # overlay filters (scaled)
    for b,(fw,ft) in FILTERS.items():
        ft_norm = ft / np.max(ft)
        plt.plot(fw.value, ft_norm * np.max(f.value)*0.7, alpha=0.35)

    plt.xlabel("Wavelength (Å)")
    plt.ylabel("fλ (erg/s/cm²/Å)")
    plt.tight_layout()
    name = f"sed_{uuid.uuid4().hex}.png"
    plt.savefig(os.path.join(STATIC_DIR, name), dpi=140)
    plt.close()
    return name

# -------------------------------------------------------------
# PLOT RESIDUALS
# -------------------------------------------------------------
def plot_residuals(comp):
    bands = ["u","g","r","i","z"]
    vals=[]
    for b in bands:
        d = comp[b]["diff"]
        vals.append(0 if d=="N/A" else float(d))

    plt.figure(figsize=(5,3))
    plt.bar(bands, vals)
    plt.axhline(0,color="black")
    plt.ylabel("Synthetic - SDSS (mag)")
    plt.tight_layout()

    name = f"res_{uuid.uuid4().hex}.png"
    plt.savefig(os.path.join(STATIC_DIR,name), dpi=140)
    plt.close()
    return name

# -------------------------------------------------------------
# MAIN ROUTES
# -------------------------------------------------------------
@app.route("/")
def index():
    return render_template("index.html")

@app.route("/process", methods=["POST"])
def process():

    clean_static_folder()
    input_type = request.form.get("input_type")

    if not input_type:
        return render_template("index.html", error="Please select an input type.")

    result=None
    sed_plot=None
    comparison=None
    res_plot=None
    total_flux_val=None

    # =========================================================
    # LUMINOSITY MODE
    # =========================================================
    if input_type=="lum":
        L=float(request.form.get("luminosity") or 0)
        D=float(request.form.get("distance") or 0)
        T=float(request.form.get("temperature") or 0)

        if L<=0 or D<=0:
            return render_template("index.html", error="Enter valid luminosity and distance.")

        L_cgs=(L*u.W).to(u.erg/u.s).value
        D_cm=(D*u.m).to(u.cm).value
        Ftot=L_cgs/(4*np.pi*D_cm**2)
        total_flux_val=Ftot

        if T>0:
            w,B=blackbody_shape(T)
            scale=Ftot/np.trapz(B,w.value)
            f=B*scale*(u.erg/(u.s*u.cm**2*u.AA))
            w_qty=w
        else:
            w_vals=np.linspace(3000,11000,2000)
            w_qty=w_vals*u.AA
            f=np.ones_like(w_vals)*(Ftot/(11000-3000))*(u.erg/(u.s*u.cm**2*u.AA))

    # =========================================================
    # FLUX MODE
    # =========================================================
    elif input_type=="flux":
        F=float(request.form.get("flux_value") or 0)
        if F<=0:
            return render_template("index.html", error="Enter valid flux.")
        w_vals=np.linspace(3000,11000,2000)
        w_qty=w_vals*u.AA
        f=np.ones_like(w_vals)*(F/(11000-3000))*(u.erg/(u.s*u.cm**2*u.AA))
        total_flux_val=F

    # =========================================================
    # MAGNITUDE MODE
    # =========================================================
    elif input_type=="mag":
        m=float(request.form.get("magnitude") or 0)
        if m==0:
            return render_template("index.html", error="Enter valid magnitude.")
        f_nu=10**(-0.4*(m+48.60))*(u.erg/(u.s*u.cm**2*u.Hz))
        w_vals=np.linspace(3000,11000,2000)*u.AA
        f=(f_nu * c.to(u.AA/u.s) / (w_vals**2)).to(u.erg/(u.s*u.cm**2*u.AA))
        w_qty=w_vals
        total_flux_val=np.trapz(f.value, w_vals.value)

    # =========================================================
    # SED (FITS OR TEXT) MODE — FIXED VERSION
    # =========================================================
    elif input_type=="sed":

        url=request.form.get("fits_url","").strip()
        sedfile=request.files.get("sedfile")

        # ---- SOURCE ----
        if url:
            try:
                data=download_fits_url(url)
                lam, flux_raw=read_sdss_fits(data)
            except:
                return render_template("index.html", error="Invalid FITS URL.")
        elif sedfile and sedfile.filename:
            name=sedfile.filename.lower()
            if name.endswith(".fits") or name.endswith(".fit"):
                lam, flux_raw=read_sdss_fits(io.BytesIO(sedfile.read()))
            else:
                lam, flux_raw=read_text_sed(io.BytesIO(sedfile.read()))
                w_qty,f=lam,flux_raw
        else:
            return render_template("index.html", error="Provide SED file or FITS URL.")

        # ---- SDSS UNIT FIX ----
        f = flux_raw * 1e-17      # Convert to erg/s/cm²/Å per steradian
        f = f * (4*np.pi)         # Convert per steradian → per cm²

        # ---- FIBER APERTURE CORRECTION ----
        fiber_area = np.pi*(1.5**2)  # 3" fiber
        f = f / fiber_area

        # ---- OPTIONAL r-BAND NORMALIZATION ----
        sdss_r = request.form.get("sdss_r","")
        if sdss_r not in ["", None]:
            try:
                sdss_r_val=float(sdss_r)
                syn_r=synthetic_abmag(lam, f*(u.erg/(u.s*u.cm**2*u.AA)), r_w, r_T)
                if syn_r is not None:
                    scale=10**(-0.4*(sdss_r_val - syn_r))
                    f = f * scale
            except:
                pass

        w_qty=lam
        f = f*(u.erg/(u.s*u.cm**2*u.AA))
        total_flux_val = np.trapz(f.value, w_qty.value)

    # =========================================================
    # COMPUTE SYNTHETIC MAGNITUDES
    # =========================================================
    mu = synthetic_abmag(w_qty, f, u_w, u_T)
    mg = synthetic_abmag(w_qty, f, g_w, g_T)
    mr = synthetic_abmag(w_qty, f, r_w, r_T)
    mi = synthetic_abmag(w_qty, f, i_w, i_T)
    mz = synthetic_abmag(w_qty, f, z_w, z_T)

    def fmt(x): return "N/A" if x is None else f"{x:.3f}"

    result = {
        "F": f"{total_flux_val:.3e}",
        "mags": f"u={fmt(mu)}, g={fmt(mg)}, r={fmt(mr)}, i={fmt(mi)}, z={fmt(mz)}",
        "colors": {
            "u-g": None if (mu is None or mg is None) else round(mu-mg,3),
            "g-r": None if (mg is None or mr is None) else round(mg-mr,3),
            "r-i": None if (mr is None or mi is None) else round(mr-mi,3),
            "i-z": None if (mi is None or mz is None) else round(mi-mz,3),
        }
    }

    # SED PLOT
    sed_plot = plot_sed(w_qty, f)

    # =========================================================
    # SDSS COMPARISON
    # =========================================================
    sdss = {
        "u":request.form.get("sdss_u",""),
        "g":request.form.get("sdss_g",""),
        "r":request.form.get("sdss_r",""),
        "i":request.form.get("sdss_i",""),
        "z":request.form.get("sdss_z","")
    }

    if any(sdss.values()):
        comparison={}
        syn={"u":mu,"g":mg,"r":mr,"i":mi,"z":mz}

        for b in ["u","g","r","i","z"]:
            try:
                catalog = float(sdss[b]) if sdss[b] else None
            except:
                catalog = None

            s = syn[b]

            comparison[b] = {
                "syn": fmt(s),
                "cat": "N/A" if catalog is None else f"{catalog:.3f}",
                "diff": "N/A" if (catalog is None or s is None) else f"{(s-catalog):+.3f}"
            }

        res_plot = plot_residuals(comparison)

    return render_template("index.html", result=result,
                           sed_plot=sed_plot,
                           comparison=comparison,
                           res_plot=res_plot)

# -------------------------------------------------------------
# MAIN
# -------------------------------------------------------------
if __name__=="__main__":
    app.run(debug=True)
