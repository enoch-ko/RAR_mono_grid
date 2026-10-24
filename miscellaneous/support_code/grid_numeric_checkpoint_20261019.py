"""Numeric I/O checkpoints only; no scientific evaluation or automatic resume.
Partial streams are never final results. Completion receipts are published by
exclusive hard-link after fsync; .pending/.partial aliases are retained, not overwritten.
"""
import csv, hashlib, json, math, os
from pathlib import Path

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def json_value(value):
    if isinstance(value,dict):return {str(k):json_value(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [json_value(v) for v in value]
    if hasattr(value,"tolist"):return json_value(value.tolist())
    if hasattr(value,"item"):return json_value(value.item())
    if isinstance(value,float) and not math.isfinite(value):return None
    return value

def publish_json(path,value):
    path=Path(path)
    if path.exists():raise FileExistsError(path)
    pending=path.with_name(path.name+".pending")
    with pending.open("x") as stream:
        json.dump(json_value(value),stream,indent=2,allow_nan=False)
        stream.write("\n");stream.flush();os.fsync(stream.fileno())
    os.link(pending,path)  # Atomic exclusive publication; retain our staging alias.
    sync_directory(path.parent)

def sync_directory(path):
    fd=os.open(path,os.O_RDONLY)
    try:os.fsync(fd)
    finally:os.close(fd)

class GridCheckpoint:
    def __init__(self,path,identity,axes):
        self.path=Path(path)
        self.path.mkdir(parents=True,exist_ok=False)
        self.identity=json_value(identity)
        self.axes=json_value(axes)
        self.seen=set();self.streams={};self.writers={};self.header=None
        self.published={};self.grids_ready=False;self.full_ready=False
        publish_json(self.path/"context.json",dict(schema="grid-numeric-checkpoint-v1",
            status="partial",identity=self.identity,axes=self.axes,
            checkpoint_policy="completed conditional+nuisance cells in original evaluation order; no automatic resume",
            nan_policy="CSV preserves original nan/inf tokens and validity/reasons; JSON nonfinite values become null"))
    def _open_streams(self,fields):
        self.header=list(fields)
        for kind in ("conditional","marginalised_gas_bulge","profile_grid"):
            stream=(self.path/(kind+".partial.csv")).open("x",newline="")
            writer=csv.DictWriter(stream,fieldnames=self.header)
            writer.writeheader();stream.flush()
            self.streams[kind]=stream;self.writers[kind]=writer
    def append_cell(self,rows):
        if self.grids_ready:raise RuntimeError("checkpoint grid already published")
        kinds={row["probability_kind"] for row in rows}
        if len(rows)!=2 or kinds!={"conditional","marginalised_gas_bulge"}:
            raise ValueError("exact conditional+nuisance pair required")
        keys={(int(row["grid_i"]),int(row["grid_j"])) for row in rows}
        if len(keys)!=1:raise ValueError("cell pair coordinates differ")
        i,j=next(iter(keys))
        if (i,j) in self.seen:raise ValueError("duplicate checkpoint cell")
        if not 0<=i<len(self.axes["slopes"]) or not 0<=j<len(self.axes["log10_mean_ml"]):
            raise ValueError("cell outside axes")
        for row in rows:
            if row["galaxy"]!=self.identity["galaxy"] or row["profile"]!="log_linear":
                raise ValueError("checkpoint galaxy/profile mismatch")
            if (float(row["log10_upsilon_disk"])!=self.axes["log10_mean_ml"][j] or
                    float(row["upsilon_disk_mean"])!=self.axes["mean_ml"][j] or
                    float(row["disk_slope"])!=self.axes["slopes"][i]):
                raise ValueError("checkpoint axes differ from original rows")
        if self.header is None:self._open_streams(rows[0])
        if any(list(row)!=self.header for row in rows):raise ValueError("CSV schema changed")
        for row in rows:
            self.writers[row["probability_kind"]].writerow(row)
        self.writers["profile_grid"].writerows(rows)
        # Flush each completed cell, rather than wait for a galaxy/shard to end.
        # No per-cell fsync claim: power failure or interrupted current write can lose the tail.
        for stream in self.streams.values():stream.flush()
        self.seen.add((i,j))
    def finish_grids(self,negative_grid):
        expected={(i,j) for i in range(len(self.axes["slopes"])) for j in range(len(self.axes["log10_mean_ml"]))}
        if self.seen!=expected:raise ValueError("incomplete cell coverage")
        for stream in self.streams.values():
            stream.flush();os.fsync(stream.fileno());stream.close()
        for kind in self.streams:
            pending=self.path/(kind+".partial.csv")
            final=self.path/(kind+".csv")
            os.link(pending,final)
            self.published[final.name]=sha(final)
        publish_json(self.path/"grids_complete.json",dict(status="grids_complete_not_fully_marginalised",
            identity=self.identity,context_sha256=sha(self.path/"context.json"),
            cells=len(self.seen),files=self.published,negative_profile_mask=negative_grid,
            row_schema=self.header,science_validated=False))
        self.published["grids_complete.json"]=sha(self.path/"grids_complete.json")
        self.grids_ready=True
    def save_fully(self,row):
        if not self.grids_ready or self.full_ready or row is None:raise ValueError("grid/full stage mismatch")
        final=self.path/"fully_marginalised.csv"
        pending=self.path/"fully_marginalised.pending.csv"
        with pending.open("x",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(row));writer.writeheader();writer.writerow(row)
            stream.flush();os.fsync(stream.fileno())
        os.link(pending,final);sync_directory(self.path)
        self.published[final.name]=sha(final)
        self.full_ready=True
    def finish_galaxy(self,radius_coverage):
        if not self.full_ready:raise ValueError("no fully marginalised checkpoint")
        publish_json(self.path/"galaxy_complete.json",dict(status="galaxy_calculation_complete",
            identity=self.identity,context_sha256=sha(self.path/"context.json"),files=self.published,
            cells=len(self.seen),radius_coverage=radius_coverage,science_validated=False,
            resume_implemented=False))
    def fail(self,error):
        # Best effort only: quota/full disk can prevent the failure receipt itself.
        failures=[]
        for stream in self.streams.values():
            if not stream.closed:
                try:stream.flush();os.fsync(stream.fileno());stream.close()
                except OSError as e:failures.append(repr(e))
        publish_json(self.path/"failure.json",dict(status="failed_partial_checkpoint",
            identity=self.identity,cells_flushed=len(self.seen),grids_ready=self.grids_ready,
            fully_marginalised_saved=self.full_ready,error=repr(error),io_failures=failures))

def verify_complete(path,identity):
    """Read-only complete-product audit, NOT a resume/restart implementation."""
    path=Path(path)
    complete=json.loads((path/"galaxy_complete.json").read_text())
    context=json.loads((path/"context.json").read_text())
    if (complete.get("status")!="galaxy_calculation_complete" or
            complete["identity"]!=json_value(identity) or context["identity"]!=json_value(identity) or
            sha(path/"context.json")!=complete["context_sha256"]):
        raise ValueError("checkpoint identity/context mismatch")
    required={"conditional.csv","marginalised_gas_bulge.csv","profile_grid.csv",
        "fully_marginalised.csv","grids_complete.json"}
    if set(complete["files"])!=required:raise ValueError("incomplete checkpoint file list")
    for name,digest in complete["files"].items():
        if Path(name).name!=name or sha(path/name)!=digest:raise ValueError("checkpoint hash mismatch")
    axes=context["axes"]
    expected={(str(i),str(j)) for i in range(len(axes["slopes"])) for j in range(len(axes["log10_mean_ml"]))}
    for name,kind in (("conditional.csv","conditional"),("marginalised_gas_bulge.csv","marginalised_gas_bulge")):
        with (path/name).open() as stream:rows=list(csv.DictReader(stream))
        if len(rows)!=len(expected) or {(row["grid_i"],row["grid_j"]) for row in rows}!=expected:
            raise ValueError("checkpoint cell coverage differs")
        if any(row["probability_kind"]!=kind or row["galaxy"]!=identity["galaxy"] for row in rows):
            raise ValueError("checkpoint row identity differs")
    return complete
