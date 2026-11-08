from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from project import stage
root = stage('score')
input_path = root / "inputs" / "all_scans.csv"
train_path = root / "inputs" / "classifier_input.csv"
params_path = root / "results" / "M3_parameters.csv"
fold_path = root / "protocol" / "outer_patient_folds.csv"

full = pd.read_csv(input_path, encoding="utf-8-sig", low_memory=False)
train = pd.read_csv(train_path, encoding="utf-8-sig", low_memory=False)
params = pd.read_csv(params_path, encoding="utf-8-sig", low_memory=False)
folds = pd.read_csv(fold_path, encoding="utf-8-sig", low_memory=False, dtype={"PTID": str})
full["PTID"] = full["PTID"].astype(str)
train["PTID"] = train["PTID"].astype(str)

jsd_cols = [f"JSD_{i:04d}" for i in range(1, 1892)]
assert not full.duplicated(["PTID","MRI_IMAGE_ID"]).any()
assert full["PTID"].nunique() > 1
assert len(train) == train["PTID"].nunique()
assert len(jsd_cols) == 1891 and all(c in full.columns for c in jsd_cols) and all(c in train.columns for c in jsd_cols)
assert set(train["PTID"]).issubset(set(full["PTID"]))
assert set(params["MODEL"]) == {"M3"}
assert params.shape[0] == 25 * (1891 + 6)
assert np.isfinite(full[jsd_cols].to_numpy(dtype=float)).all()
assert np.isfinite(train[jsd_cols].to_numpy(dtype=float)).all()

train = train.set_index("PTID", drop=False)
train_ids = set(train.index)
fold_map = {(str(r.PTID), int(r.REPEAT)): int(r.OUTER_FOLD) for r in folds.itertuples(index=False)}
assert len(fold_map) == len(train) * 5
for rep in range(1, 6):
    assert set(folds.loc[folds.REPEAT.eq(rep), "PTID"].astype(str)) == train_ids

full_ids = full["PTID"].to_numpy(dtype=str)
unique_ids = pd.Index(pd.unique(full_ids))
assignments = {}
for ptid in unique_ids:
    for rep in range(1, 6):
        if ptid in train_ids:
            fold = fold_map[(ptid, rep)]
            source = "CNAD_HELD_OUT"
        else:
            token = f"2026102403|{rep}|{ptid}".encode()
            fold = 1 + int(hashlib.sha256(token).hexdigest()[:16], 16) % 5
            source = "NEVER_TRAINED_HASH_FOLD"
        assignments[(ptid, rep)] = (fold, source)

X = full[jsd_cols].to_numpy(dtype=np.float64)
X_train = train[jsd_cols].to_numpy(dtype=np.float64)
rows = []
refs = []
all_raw = np.zeros((len(full), 5), dtype=np.float64)
all_cnref = np.zeros((len(full), 5), dtype=np.float64)
all_fold = np.zeros((len(full), 5), dtype=np.int16)
all_source = np.empty((len(full), 5), dtype=object)
for rep in range(1, 6):
    for fold in range(1, 6):
        p = params[(params.REPEAT == rep) & (params.OUTER_FOLD == fold)].copy()
        assert len(p) == 1897
        p = p.set_index("FEATURE")
        assert p.index[0] == "(Intercept)"
        pp = p.loc[jsd_cols]
        beta = pp["BETA"].to_numpy(dtype=np.float64)
        mean = pp["MEAN"].to_numpy(dtype=np.float64)
        scale = pp["SCALE"].to_numpy(dtype=np.float64)
        assert np.all(np.isfinite(beta)) and np.all(np.isfinite(mean)) and np.all(scale > 0)
        raw_full = ((X - mean) / scale) @ beta
        raw_train = ((X_train - mean) / scale) @ beta
        train_cn = train["Y"].to_numpy(dtype=int) == 0
        train_fold = folds.loc[(folds.REPEAT == rep) & (folds.OUTER_FOLD != fold), "PTID"].astype(str).to_numpy()
        train_cn_ids = [v for v in train_fold if int(train.loc[v, "Y"]) == 0]
        cn_ix = train.index.get_indexer(train_cn_ids)
        assert np.all(cn_ix >= 0) and len(cn_ix) > 0
        cn_mean = float(raw_train[cn_ix].mean())
        cn_sd = float(raw_train[cn_ix].std(ddof=0))
        assert np.isfinite(cn_mean) and np.isfinite(cn_sd) and cn_sd > 0
        for i, ptid in enumerate(full_ids):
            fold_i, source_i = assignments[(ptid, rep)]
            if fold_i == fold:
                all_raw[i, rep - 1] = raw_full[i]
                all_cnref[i, rep - 1] = (raw_full[i] - cn_mean) / cn_sd
                all_fold[i, rep - 1] = fold
                all_source[i, rep - 1] = source_i
        n_applied = sum(1 for ptid in unique_ids if assignments[(ptid, rep)][0] == fold)
        overlap = len(set([ptid for ptid in unique_ids if assignments[(ptid, rep)][0] == fold]) & set(train_fold))
        assert overlap == 0
        refs.append({"REPEAT": rep, "OUTER_FOLD": fold, "N_TRAIN": len(train_fold), "N_TRAIN_CN": len(train_cn_ids), "CN_RAW_MEAN": cn_mean, "CN_RAW_POP_SD": cn_sd, "N_APPLIED_SUBJECTS": n_applied, "N_APPLIED_MRI": int(sum(1 for ptid in full_ids if assignments[(ptid, rep)][0] == fold)), "N_TRAIN_APPLIED_OVERLAP": overlap, "JSD_FEATURES": 1891, "JSD_BETA_L2": float(np.linalg.norm(beta))})

for i, ptid in enumerate(full_ids):
    for rep in range(5):
        fold_i, source_i = assignments[(ptid, rep + 1)]
        assert all_fold[i, rep] == fold_i and all_source[i, rep] == source_i
        assert np.isfinite(all_raw[i, rep]) and np.isfinite(all_cnref[i, rep])

score_cnref = all_cnref.mean(axis=1)
first_idx = full.sort_values(["PTID", "MRI_DATE", "MRI_IMAGE_ID"], kind="mergesort").drop_duplicates("PTID", keep="first").index.to_numpy()
assert len(first_idx) == full["PTID"].nunique()
baseline_mean = float(score_cnref[first_idx].mean())
baseline_sd = float(score_cnref[first_idx].std(ddof=1))
assert np.isfinite(baseline_mean) and np.isfinite(baseline_sd) and baseline_sd > 0
score_z = (score_cnref - baseline_mean) / baseline_sd

out = full.copy()
out["JSD_SCORE_CNREF"] = score_cnref
out["JSD_SCORE_Z"] = score_z
for rep in range(1, 6):
    out[f"JSD_SCORE_R{rep}"] = all_cnref[:, rep - 1]
    out[f"JSD_RAW_R{rep}"] = all_raw[:, rep - 1]
    out[f"JSD_APPLIED_FOLD_R{rep}"] = all_fold[:, rep - 1]
    out[f"JSD_APPLIED_SOURCE_R{rep}"] = all_source[:, rep - 1]
out["JSD_BASELINE_ROW"] = False
out.loc[first_idx, "JSD_BASELINE_ROW"] = True

result_dir = root / "results"
audit_dir = root / "audit"
result_dir.mkdir(parents=True, exist_ok=True)
audit_dir.mkdir(parents=True, exist_ok=True)
out.to_csv(result_dir / "scored_MRI.csv", index=False, encoding="utf-8-sig", float_format="%.17g")
out.loc[first_idx, ["PTID", "MRI_IMAGE_ID", "MRI_DATE", "DX_STD", "JSD_SCORE_CNREF", "JSD_SCORE_Z"]].sort_values("PTID").to_csv(result_dir / "baseline_scores.csv", index=False, encoding="utf-8-sig", float_format="%.17g")
pd.DataFrame(refs).to_csv(audit_dir / "score_model_summary.csv", index=False, encoding="utf-8-sig", float_format="%.17g")
application_rows = []
for ptid in unique_ids:
    for rep in range(1, 6):
        fold_i, source_i = assignments[(ptid, rep)]
        n_mri = int((full_ids == ptid).sum())
        application_rows.append({"PTID": ptid, "REPEAT": rep, "OUTER_FOLD": fold_i, "SOURCE": source_i, "N_MRI": n_mri, "IN_SOURCE_TRAINING": ptid in train_ids})
pd.DataFrame(application_rows).to_csv(audit_dir / "application_assignments.csv", index=False, encoding="utf-8-sig")

summary = {
    "n_mri": int(len(out)),
    "n_participants": int(out.PTID.nunique()),
    "n_baseline": int(len(first_idx)),
    "n_jsd_features": 1891,
    "n_models": 25,
    "source_training_participants": int(len(train_ids)),
    "source_never_trained_participants": int(len(unique_ids) - len(train_ids)),
    "baseline_mean_cnref": baseline_mean,
    "baseline_sample_sd_cnref": baseline_sd,
    "score_min": float(np.nanmin(score_z),),
    "score_max": float(np.nanmax(score_z),),
    "source_model": "M3",
    "scoring": "five repeat-specific held-out outer-fold components averaged per MRI; fixed model per participant across visits",
    "training_overlap": 0,
    "all_scores_finite": bool(np.isfinite(score_z).all()),
}
(audit_dir / "score_scale.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False))
