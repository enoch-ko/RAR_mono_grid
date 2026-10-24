"""Named, header-validated SPARC rotmod columns; velocities in km/s, radii in kpc.
The published Vgas template already includes helium. Do not scale it twice.
"""
from pathlib import Path
import hashlib
import numpy as np

ROTMOD_COLUMNS=("Rad","Vobs","errV","Vgas","Vdisk","Vbul","SBdisk","SBbul")
ROTMOD_INDEX={name:i for i,name in enumerate(ROTMOD_COLUMNS)}
ROTMOD_PARSER_ID="sparc_rotmod_header_v1_gas3_disk4"


def parser_provenance():
    return {"id":ROTMOD_PARSER_ID,"columns":ROTMOD_INDEX.copy(),
            "source_file":Path(__file__).name,
            "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def read_rotmod(path):
    path=Path(path)
    header=next((line.lstrip("# ").split() for line in path.read_text().splitlines()
                 if line.startswith("#") and "Vobs" in line),None)
    if header!=list(ROTMOD_COLUMNS):
        raise ValueError("rotmod header must declare Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul")
    values=np.loadtxt(path,ndmin=2)
    if values.ndim!=2 or values.shape[1]!=len(ROTMOD_COLUMNS):
        raise ValueError("rotmod requires exactly eight declared columns")
    return values


def read_column(values,name,mask=None):
    values=np.asarray(values)
    if values.ndim!=2 or values.shape[1]!=len(ROTMOD_COLUMNS):
        raise ValueError("rotmod array requires eight columns")
    if name not in ROTMOD_INDEX:raise ValueError("unknown rotmod column: "+str(name))
    result=values[:,ROTMOD_INDEX[name]]
    return result if mask is None else result[mask]


def disc_reference_velocity(values,mask,disc_ml):
    if not np.isfinite(disc_ml) or disc_ml<=0:raise ValueError("positive finite disc M/L required")
    return read_column(values,"Vdisk",mask)*np.sqrt(disc_ml)
