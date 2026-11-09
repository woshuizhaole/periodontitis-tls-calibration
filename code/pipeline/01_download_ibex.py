#!/usr/bin/env python3
import json, os, time, hashlib
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests

ROOT="/data/workplace/wyx/perio_tls_20261107"
OUT=os.path.join(ROOT,"data/raw/ibex")
LOG=os.path.join(ROOT,"logs/downloads/ibex_download.log")
os.makedirs(OUT,exist_ok=True); os.makedirs(os.path.dirname(LOG),exist_ok=True)
DOI="doi:10.7910/DVN/3PFU0D"; API="https://dataverse.harvard.edu/api"
H={"User-Agent":"Mozilla/5.0"}

def log(m):
    line=f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"; print(line,flush=True)
    open(LOG,"a").write(line+"\n")

def main():
    r=requests.get(f"{API}/datasets/:persistentId/?persistentId={DOI}",headers=H,timeout=120); r.raise_for_status()
    files=r.json()["data"]["latestVersion"]["files"]
    log(f"files={len(files)}")
    jobs=[]
    for f in files:
        df=f["dataFile"]; name=df["filename"]; size=df.get("filesize",0)
        if not name.endswith(".ims"): continue
        jobs.append((df["id"],name,size))
    log(f"ims jobs={len(jobs)} total={sum(j[2] for j in jobs)/1e9:.2f} GB")
    def dl(job):
        fid,name,size=job; out=os.path.join(OUT,name)
        if os.path.exists(out) and os.path.getsize(out)==size: return name,size,"exists"
        for attempt in range(1,6):
            try:
                with requests.get(f"{API}/access/datafile/{fid}",headers=H,stream=True,timeout=600) as r:
                    r.raise_for_status(); tmp=out+".part"
                    with open(tmp,"wb") as fh:
                        for chunk in r.iter_content(1024*1024*8):
                            if chunk: fh.write(chunk)
                    os.replace(tmp,out)
                if os.path.getsize(out)!=size: raise IOError("size mismatch")
                return name,size,"downloaded"
            except Exception as e:
                log(f"RETRY {name} {attempt}: {e}"); time.sleep(5*attempt)
        raise RuntimeError(name)
    rows=[]
    with ThreadPoolExecutor(max_workers=3) as ex:
        futs={ex.submit(dl,j):j for j in jobs}
        for i,fut in enumerate(as_completed(futs),1):
            name,size,status=fut.result(); rows.append((name,size,status)); log(f"{i}/{len(jobs)} {status} {name}")
    with open(os.path.join(ROOT,"data/manifests/ibex_files.tsv"),"w") as fh:
        fh.write("file\tsize\tstatus\n")
        for row in rows: fh.write("\t".join(map(str,row))+"\n")
    log("done")

if __name__=="__main__":
    main()
