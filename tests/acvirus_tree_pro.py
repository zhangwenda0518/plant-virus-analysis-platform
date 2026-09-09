#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ACVirus Tree-Pro (Streamlined Edition) — 按科建树 + 属级分色 + 共线性双面板"""
from __future__ import annotations
import argparse
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
from matplotlib.colors import LinearSegmentedColormap, Normalize

COMMAND_LOG_NAME = "tree_pro_commands.jsonl"
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})

def find_executable(names): 
    for name in names:
        path = shutil.which(name)
        if path: return path
    return None

def run_cmd(cmd, label, log_dir, stdout_file=None):
    cmd = [str(x) for x in cmd]
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / COMMAND_LOG_NAME
    start = datetime.now().isoformat(timespec="seconds")
    stdout_handle = stdout_file.open("w") if stdout_file else subprocess.PIPE
    print(f"[RUN] {label} -> {' '.join(cmd[:6])}{' ...' if len(cmd) > 6 else ''}")
    res = subprocess.run(cmd, stdout=stdout_handle, stderr=subprocess.PIPE, text=True)
    if stdout_file: stdout_handle.close()
    with log_file.open("a", encoding="utf-8") as h:
        h.write(json.dumps({"label": label, "cmd": shlex.join(cmd), "start": start,
                            "returncode": res.returncode, "stderr": res.stderr}) + "\n")
    if res.returncode != 0:
        print(f"[ERROR] {label} failed:\n{res.stderr}", file=sys.stderr)
        raise subprocess.CalledProcessError(res.returncode, cmd, stderr=res.stderr)
    return res

def get_hierarchical_accessions(df_taxa, target_rank, target_name, mode,
                                target_ratio=1.0, other_count=3):
    """同科内智能抽样 (macro/genus/lineage), 防止全量抓取"""
    df = df_taxa.copy()
    for col in ['Realm','Order','Class','Family','Genus','Species']:
        if col in df.columns:
            df[col] = df[col].fillna("").astype(str).str.strip()
    sampled_accs = set()
    target_info = df[df[target_rank].str.lower() == target_name.lower()]
    if target_info.empty:
        raise ValueError(f"Cannot find {target_rank} '{target_name}' in database taxonomy!")
    t_family = target_info.iloc[0].get('Family','')
    t_genus  = target_info.iloc[0].get('Genus','')
    print(f"\n[SAMPLING] Mode: {mode.upper()} | Target: {target_rank}={target_name} (Family: {t_family})")
    if mode == "macro":
        fam_df = df[df['Family'].str.lower()==t_family.lower()]
        for genus, group in fam_df.groupby('Genus'):
            if genus == "" or pd.isna(genus): continue
            if str(genus).lower()==str(t_genus).lower():
                n = max(1, int(len(group)*target_ratio))
                sampled = group.sample(n=n, random_state=42)
            else:
                sampled = group.sample(n=min(len(group), other_count), random_state=42)
            sampled_accs.update(sampled['Virus GENBANK accession'].tolist())
    elif mode == "genus":
        target_group = df[df['Genus'].str.lower()==str(t_genus).lower()]
        sampled_accs.update(target_group['Virus GENBANK accession'].tolist())
        fam_other = df[(df['Family'].str.lower()==t_family.lower()) & (df['Genus'].str.lower()!=str(t_genus).lower())]
        if not fam_other.empty:
            for _, group in fam_other.groupby('Genus'):
                sampled_accs.update(group.sample(n=min(len(group), other_count), random_state=42)['Virus GENBANK accession'].tolist())
    elif mode == "lineage":
        sp_name = target_name if target_rank=='Species' else target_info.iloc[0].get('Species','')
        target_group = df[df['Species'].str.lower()==str(sp_name).lower()]
        sampled_accs.update(target_group['Virus GENBANK accession'].tolist())
        same_genus_other = df[(df['Genus'].str.lower()==str(t_genus).lower()) & (df['Species'].str.lower()!=str(sp_name).lower())]
        if not same_genus_other.empty:
            for _, group in same_genus_other.groupby('Species'):
                sampled_accs.update(group.sample(n=1, random_state=42)['Virus GENBANK accession'].tolist())
    final_accs = [str(acc).split(';')[0].strip() for acc in sampled_accs if pd.notna(acc) and str(acc).strip()]
    df['Clean_Acc'] = df['Virus GENBANK accession'].apply(lambda x: str(x).split('.')[0].split(';')[0].strip() if pd.notna(x) else '')
    taxa_meta = df.set_index('Clean_Acc')[['Species','Genus','Family','Order']].to_dict('index')
    print(f"  -> Sampled {len(final_accs)} representative genomes from database.")
    return final_accs, taxa_meta

def parse_best_model(model_out_file):
    if not Path(model_out_file).exists(): return None
    with open(model_out_file,'r') as f:
        for line in f:
            if line.strip().startswith("Best model according to"):
                for next_line in f:
                    if next_line.strip().startswith("Model:"):
                        return next_line.strip().split()[1]
    return None

def run_phylogeny_engine(aln_file, outdir, threads, engine_preference="iqtree", bootstrap=1000,
                          rooting="unrest", aln_fasta=None):
    """建树引擎: (默认) IQ-TREE + UNREST 非可逆模型 → 天然有根树 (无外群定根)
       也可用 MFP (ModelFinder 可逆模型, 无根) + 后续定根."""
    tree_file = outdir / "phylogeny.treefile"
    modeltest = find_executable(["modeltest-ng"]); raxml = find_executable(["raxml-ng"])
    iqtree = find_executable(["iqtree2","iqtree","iqtree3"])
    fasttree = find_executable(["FastTree","FastTreeDbl","fasttree"])
    # FastTree 优先: 秒级出树, SH-like 支持值 (0-100) 与 bootstrap 同尺度
    if engine_preference == "fasttree" or (engine_preference == "auto" and fasttree):
        if not fasttree: raise FileNotFoundError("--engine fasttree requested but not found!")
        # 喂 FASTA: trimAl 的交错 phylip 后续块顶格无名字区, FastTree 会误判
        ft_in = aln_fasta if (aln_fasta and Path(aln_fasta).exists()) else aln_file
        print(f"[TREE] FastTree GTR+Gamma ({Path(ft_in).name}, SH-like supports, unrooted)")
        run_cmd([fasttree,"-nt","-gtr","-gamma","-quiet",str(ft_in)],"FastTree",outdir,
                stdout_file=tree_file)
        return tree_file
    # 大比对 (>80 序列) 跳过 ModelTest+RAxML-NG (bootstrap 极慢易超时), 直接 IQ-TREE
    try:
        with open(aln_file) as _f: n_seq_aln = int(_f.readline().split()[0])
    except Exception:
        n_seq_aln = 0
    use_raxml = (engine_preference == "raxml") or (
        engine_preference == "auto" and bool(modeltest and raxml) and n_seq_aln <= 80)
    if use_raxml:
        print("[TREE] ModelTest-NG -> RAxML-NG (注意: RAxML 产无根树)")
        try:
            mt = outdir/"modeltest"
            run_cmd([modeltest,"-i",str(aln_file),"-d","nt","-p",str(threads),"-o",str(mt)],"ModelTest-NG",outdir)
            best_model = parse_best_model(str(mt)+".out") or "GTR+G"
            print(f"  -> Model: {best_model}")
            rp = outdir/"raxml"
            run_cmd([raxml,"--all","--msa",str(aln_file),"--model",best_model,"--threads",str(threads),
                     "--bs-trees",str(bootstrap),"--prefix",str(rp),"--redo"],"RAxML-NG",outdir)
            bt = Path(str(rp)+".raxml.support")
            if not bt.exists(): bt = Path(str(rp)+".raxml.bestTree")
            if bt.exists():
                shutil.copy2(bt, tree_file); return tree_file
        except Exception as e:
            print(f"[WARNING] RAxML failed ({e}), IQ-TREE fallback.")
    if not iqtree: raise FileNotFoundError("No RAxML/IQ-TREE!")
    # 默认: UNREST 非可逆模型 → 天然有根树 (无外群定根)
    if rooting == "unrest":
        print("[TREE] IQ-TREE UNREST (non-reversible) -> rooted tree (no outgroup needed)")
        ip = outdir/"iqtree_run_unrest"
        run_cmd([iqtree,"-s",str(aln_file),"-m","UNREST","--model-joint","UNREST","-B",str(bootstrap),
                 "-T",str(threads),"--prefix",str(ip),"-redo"],"IQ-TREE UNREST",outdir)
        bt = Path(str(ip)+".treefile")
    else:
        print("[TREE] IQ-TREE ModelFinder (unrooted, MFP)")
        ip = outdir/"iqtree_run"
        run_cmd([iqtree,"-s",str(aln_file),"-m","MFP","-B",str(bootstrap),"-T",str(threads),
                 "--prefix",str(ip),"-redo"],"IQ-TREE",outdir)
        bt = Path(str(ip)+".treefile")
    shutil.copy2(bt, tree_file); return tree_file

def parse_prodigal_gff(gff_file, fasta_file):
    fasta_lens = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta_file,"fasta")}
    rows=[]; pc=defaultdict(int)
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p=line.strip().split("\t")
            if len(p)>=9 and p[2]=="CDS":
                seqid,start,end = p[0].strip(),int(p[3]),int(p[4])
                if start>end: start,end=end,start
                pc[seqid]+=1
                rows.append({"nucl_id":seqid,"protein":f"{seqid}_{pc[seqid]}","start":start,"end":end,
                             "nucl_length":fasta_lens.get(seqid,0),"strand":p[6] if p[6] in['+','-'] else '+'})
    return pd.DataFrame(rows)

def prepare_synteny_and_diamond(fasta_file, outdir, threads):
    prodigal=find_executable(["prodigal"]); diamond=find_executable(["diamond"])
    if not (prodigal and diamond): raise FileNotFoundError("prodigal/diamond required")
    faa=outdir/"all_proteins.faa"; gff=outdir/"all_proteins.gff"
    dmnd=outdir/"proteins.dmnd"; dt=outdir/"diamond_all_vs_all.tsv"
    run_cmd([prodigal,"-i",str(fasta_file),"-a",str(faa),"-f","gff","-p","meta","-o",str(gff),"-q"],"Prodigal",outdir)
    pos_df=parse_prodigal_gff(gff,fasta_file)
    run_cmd([diamond,"makedb","--in",str(faa),"-d",str(dmnd),"--threads",str(threads),"--quiet"],"DIAMOND makedb",outdir)
    run_cmd([diamond,"blastp","-q",str(faa),"-d",str(dmnd),"-p",str(threads),"-o",str(dt),
             "--evalue","1e-5","--outfmt","6","qseqid","sseqid","pident","length","bitscore","--quiet"],"DIAMOND blastp",outdir)
    sim=pd.read_csv(dt,sep="\t",header=None,names=["protein1","protein2","identity","length","bitscore"])
    sim=sim[sim["protein1"]!=sim["protein2"]].copy()
    return pos_df, sim

def _root_path_len(tree, node):
    d=0.0; cur=node
    while cur is not tree.root:
        d += (cur.branch_length or 0.0)
        cur = _parent_of(tree, cur)
    return d

def _parent_of(tree, node):
    for clade in tree.get_nonterminals():
        if node in clade.clades:
            return clade
    return tree.root

def _dist_between(tree, a, b):
    # 双节点到根路径长之和近似
    return _root_path_len(tree,a) + _root_path_len(tree,b)

def _root_at_midpoint(tree, far1=None, far2=None):
    """mid 定根: 直径(最远两片叶), 选第二个为 outgroup 重根定根"""
    try:
        terms = tree.get_terminals()
        if len(terms) < 2:
            return tree
        if far1 is None:
            far1 = max(terms, key=lambda t: _root_path_len(tree, t))
        if far2 is None:
            far2 = max(terms, key=lambda t: _dist_between(tree, far1, t))
        if far1.name != far2.name:
            tree.root_with_outgroup(far2)
        return tree
    except Exception as e:
        return tree

def draw_tree_synteny_figure(tree_file, pos_df, sim_df, taxa_meta, out_prefix,
                             ignore_branch_length=False, min_identity=30.0, show_support=70,
                             keep_tree_order=False):
    tree=Phylo.read(tree_file,"newick")
    if ignore_branch_length:
        for c in tree.find_clades(): c.branch_length=1.0
    # mid 定根: UNREST 非可逆模型已天然定根, 不再手动重根; 仅非 UNREST 引擎调用时
    # _root_at_midpoint(tree) 已移除 - UNREST 建树直接得到有根树
    try:
        if not keep_tree_order:
            tree.ladderize(reverse=True)
    except Exception:
        pass
    terminals=tree.get_terminals(); num_taxa=len(terminals)
    genome_order=[n.name for n in terminals]
    y_coords={name:i for i,name in enumerate(genome_order)}
    fig=plt.figure(figsize=(19,max(6,num_taxa*0.45)),dpi=300)
    gs=fig.add_gridspec(1,3,width_ratios=[0.09,0.40,0.51],wspace=0.05)
    ax_legend=fig.add_subplot(gs[0])
    ax_tree=fig.add_subplot(gs[1])
    # 图例列（独立轴, 仅绘图例色块, 不画树）
    def assign_pos(clade,curr_x):
        clade.x=curr_x+(clade.branch_length or 0.0)
        if clade.is_terminal():
            clade.y=y_coords[clade.name]
        else:
            for child in clade.clades: assign_pos(child,clade.x)
            clade.y=np.mean([c.y for c in clade.clades])
    assign_pos(tree.root,0.0)
    max_depth=max(n.x for n in terminals)
    # 属级背景色块
    blocks=[]; curr_g=None; start_y=None
    for n in terminals:
        base_acc=n.name.split('.')[0]
        g=str(taxa_meta.get(base_acc,{}).get('Genus','Sample_Contig'))
        if g!=curr_g:
            if curr_g is not None: blocks.append((curr_g,start_y,y_coords[n.name]-1))
            curr_g,start_y=g,y_coords[n.name]
    if curr_g is not None: blocks.append((curr_g,start_y,num_taxa-1))
    cmap_bg=plt.get_cmap('Pastel1')
    gen_colors={g:cmap_bg(i%9) for i,(g,_,_) in enumerate(blocks)}
    # 分色: 我们的样本标红色斜体, 参考按属色块
    sample_colors=gen_colors
    for g,sy,ey in blocks:
        rect=patches.Rectangle((0,sy-0.4),max_depth*1.05,(ey-sy)+0.8,
                 facecolor=gen_colors.get(g,'#eeeeee'),alpha=0.35,lw=0,zorder=0)
        ax_tree.add_patch(rect)
    def draw_clade_lines(clade):
        if not clade.is_terminal():
            ys=[c.y for c in clade.clades]
            ax_tree.plot([clade.x,clade.x],[min(ys),max(ys)],color="#2c3e50",lw=1.2)
            for c in clade.clades:
                ax_tree.plot([clade.x,c.x],[c.y,c.y],color="#2c3e50",lw=1.2)
                # 标注支持值 (bootstrap/SH-like) - 仅显示 >= show_support
                supp = getattr(c, 'confidence', None)
                if supp is not None:
                    try: supp = float(supp)
                    except Exception: supp = None
                if supp is not None:
                    if supp < 1.0: supp *= 100.0  # FastTree SH-like 0-1 → 0-100 尺度
                    if supp >= float(show_support):
                        ax_tree.text(clade.x + (c.x-clade.x)/2, c.y, str(int(round(supp))),
                                     ha='center', va='center', fontsize=6.5, color='#111111',
                                     fontweight='bold', zorder=8, bbox=dict(boxstyle='round,pad=0.12',
                                     fc='white', ec='none', alpha=0.75))
                draw_clade_lines(c)
    draw_clade_lines(tree.root)
    align_x=max_depth*1.08
    for n in terminals:
        base_acc=n.name.split('.')[0]
        sp=taxa_meta.get(base_acc,{}).get('Species','')
        label=f"{sp} ({n.name})" if sp else n.name
        ax_tree.plot([n.x,align_x],[n.y,n.y],color="#bdc3c7",linestyle=":",lw=0.8)
        is_sample = base_acc not in taxa_meta
        ax_tree.text(align_x+max_depth*0.02,n.y,label,va="center",ha="left",fontsize=8,
                     fontweight="bold" if is_sample else "normal", color="#c0392b" if is_sample else "#2c3e50")
    ax_tree.set_xlim(0,max_depth*1.6); ax_tree.set_ylim(num_taxa-0.5,-0.5)
    ax_tree.axis("off"); ax_tree.set_title("Maximum-Likelihood Phylogeny",loc="left",fontsize=11,fontweight="bold",pad=12)
    # 属图例 - 独立列 (ax_legend), 顶部色块+标签竖直排列, 完全不接触树
    y_leg = 0.92
    for g,_,_ in blocks[:12]:
        ax_legend.add_patch(patches.Rectangle((0.15,y_leg), 0.35, 0.045,
                          facecolor=gen_colors.get(g,'#ccc'), edgecolor='#555',linewidth=0.6, transform=ax_legend.transAxes))
        ax_legend.text(0.55, y_leg+0.018, g, fontsize=8, va='center', ha='left', transform=ax_legend.transAxes)
        y_leg -= 0.07
    if blocks:
        ax_legend.text(0.15, 0.97, 'Genus', fontsize=9, fontweight='bold', transform=ax_legend.transAxes, va='top')
    ax_legend.axis('off')
    ax_legend.set_ylim(-0.1, 1.2)
    # 右面板 synteny
    ax_syn=fig.add_subplot(gs[2])
    max_len=pos_df["nucl_length"].max() if not pos_df.empty else 10000
    pos_df_sub=pos_df[pos_df["nucl_id"].isin(genome_order)].copy()
    ymax=0
    for genome,y in y_coords.items():
        gl=pos_df[pos_df["nucl_id"]==genome]["nucl_length"].max()
        if pd.notna(gl):
            ymax=max(ymax,gl)
            ax_syn.hlines(y,0,gl,color="#34495e",lw=2,zorder=2)
    for row in pos_df_sub.itertuples():
        y=y_coords[row.nucl_id]
        w=max(row.end-row.start,1.0)
        rect=patches.FancyBboxPatch((row.start,y-0.15),w,0.3,boxstyle="round,pad=0.02",facecolor="#3498db",edgecolor="#2980b9",lw=0.5,zorder=4)
        ax_syn.add_patch(rect)
    prot_to_info=pos_df_sub.set_index("protein").to_dict("index")
    sim_cmap=LinearSegmentedColormap.from_list("sim",["#d5e8f7","#3498db","#2ecc71","#e67e22"])
    norm=Normalize(vmin=min_identity,vmax=100.0)
    for row in sim_df.itertuples():
        if row.identity<min_identity: continue
        p1,p2=str(row.protein1),str(row.protein2)
        if p1 in prot_to_info and p2 in prot_to_info:
            i1,i2=prot_to_info[p1],prot_to_info[p2]; g1,g2=i1["nucl_id"],i2["nucl_id"]
            if abs(y_coords[g1]-y_coords[g2])==1:
                y1,y2=y_coords[g1],y_coords[g2]
                s1=a if False else i1["start"]; e1=i1["end"]
                y1u = y1+0.15 if y1<y2 else y1-0.15
                y2u = y2-0.15 if y1<y2 else y2+0.15
                poly=patches.Polygon([(i1["start"],y1u),(i1["end"],y1u),(i2["end"],y2u),(i2["start"],y2u)],
                      closed=True,facecolor=sim_cmap(norm(row.identity)),alpha=0.35,edgecolor="none",zorder=3)
                ax_syn.add_patch(poly)
    ax_syn.set_xlim(-max_len*0.02,max_len*1.05); ax_syn.set_ylim(num_taxa-0.5,-0.5)
    ax_syn.set_xlabel("Genome Coordinates (bp)",fontsize=9,fontweight="bold")
    ax_syn.spines['top'].set_visible(False); ax_syn.spines['right'].set_visible(False); ax_syn.spines['left'].set_visible(False)
    ax_syn.set_yticks([]); ax_syn.set_title("Comparative Synteny & Homology",loc="left",fontsize=11,fontweight="bold",pad=12)
    cbar_ax=ax_syn.inset_axes([0.72,1.02,0.25,0.03])
    cbar=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=sim_cmap),cax=cbar_ax,orientation="horizontal")
    cbar.set_label("Protein Identity (%)",fontsize=7); cbar.ax.tick_params(labelsize=6)
    fig.savefig(f"{out_prefix}.png",dpi=300,bbox_inches="tight")
    fig.savefig(f"{out_prefix}.pdf",bbox_inches="tight")
    plt.close(fig)
    print(f"[PLOT] Figures: {out_prefix}.png / .pdf")

def main():
    parser=argparse.ArgumentParser(prog="acvirus_tree_pro",description="ACVirus Tree-Pro: Phylogeny & Synteny",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--mode",choices=["macro","genus","lineage","custom"],default="macro")
    parser.add_argument("--target_name",type=str,help="Target taxon name (e.g. Rhabdoviridae)")
    parser.add_argument("--target_rank",choices=["Family","Genus","Species"],default="Family")
    parser.add_argument("--contigs",type=Path,help="Input query viral contigs FASTA")
    parser.add_argument("--db_taxa",type=Path,required=True)
    parser.add_argument("--db_fasta",type=Path,required=True)
    parser.add_argument("--outdir",type=Path,required=True)
    parser.add_argument("--target_ratio",type=float,default=1.0)
    parser.add_argument("--other_count",type=int,default=3)
    parser.add_argument("--engine",choices=["auto","fasttree","raxml","iqtree"],default="auto")
    parser.add_argument("--threads",type=int,default=16)
    parser.add_argument("--bootstrap",type=int,default=1000)
    parser.add_argument("--min_identity",type=float,default=30.0)
    parser.add_argument("--ignore_branch_length",action="store_true")
    parser.add_argument("--show_support",type=float,default=70.0,help="Minimum branch support displayed (default 70)")
    parser.add_argument("--keep_tree_order",action="store_true",help="Preserve tree leaf order (no ladderizing)")
    parser.add_argument("--rooting",choices=["unrest","mfp"],default="mfp",help="Tree rooting: unrest (UNREST non-reversible, rooted) or mfp (ModelFinder, default unrooted)")
    args=parser.parse_args()
    outdir=args.outdir.resolve(); outdir.mkdir(parents=True,exist_ok=True)
    df_taxa=pd.read_csv(args.db_taxa)
    if args.mode!="custom":
        if not args.target_name: parser.error("--target_name required!")
        sampled_accs,taxa_meta=get_hierarchical_accessions(df_taxa,args.target_rank,args.target_name,args.mode,args.target_ratio,args.other_count)
    else:
        sampled_accs,taxa_meta=[],{}
    seqkit=find_executable(["seqkit"])
    if not seqkit: raise FileNotFoundError("seqkit required!")
    final_fasta=outdir/"analysis_sequences.fasta"
    acc_file=outdir/"sampled_accessions.txt"
    with acc_file.open("w") as f:
        for acc in sampled_accs: f.write(f"^{acc}(\\.[0-9]+)?\n")
    db_extracted=outdir/"db_extracted.fasta"
    if sampled_accs:
        run_cmd([seqkit,"grep","-r","-f",str(acc_file),str(args.db_fasta),"-o",str(db_extracted)],"Extract DB",outdir)
    with final_fasta.open("wb") as out_f:
        if args.contigs and args.contigs.exists():
            with args.contigs.open("rb") as inf: shutil.copyfileobj(inf,out_f)
        if db_extracted.exists():
            with db_extracted.open("rb") as inf: shutil.copyfileobj(inf,out_f)
    mafft=find_executable(["mafft"]); trimal=find_executable(["trimal"])
    aln_out=outdir/"alignment.mafft"; trimmed=outdir/"alignment.trimmed.phy"
    run_cmd([mafft,"--auto","--thread",str(args.threads),str(final_fasta)],"MAFFT",outdir,stdout_file=aln_out)
    # 统一大写 (soft-mask 小写符号会让 trimAl 报 symbol not defined)
    upper_aln=outdir/"alignment.upper.fasta"
    with aln_out.open() as _fi, upper_aln.open("w") as _fo:
        for _line in _fi: _fo.write(_line if _line.startswith(">") else _line.upper())
    aln_out=upper_aln
    run_cmd([trimal,"-in",str(aln_out),"-out",str(trimmed),"-gt","0.8","-st","0.005","-phylip"],"trimAl",outdir)
    # 保护: trim 后序列数不足则跳过 (trimal 删光全 gap 序列后不会写出文件)
    def _count_phylip(p):
        try:
            with p.open() as _f:
                return int(_f.readline().split()[0])
        except Exception:
            return 0
    n_after=_count_phylip(trimmed)
    if n_after<3:
        print(f"[SKIP] {outdir.name}: only {n_after} sequences after trimming, skip tree")
        return
    tree_file=run_phylogeny_engine(trimmed,outdir,args.threads,args.engine,args.bootstrap,rooting=getattr(args,'rooting','mfp'),aln_fasta=aln_out)
    pos_df,sim_df=prepare_synteny_and_diamond(final_fasta,outdir,args.threads)
    draw_tree_synteny_figure(tree_file,pos_df,sim_df,taxa_meta,out_prefix=outdir/"Tree_Synteny_Composite",
                             ignore_branch_length=args.ignore_branch_length,min_identity=args.min_identity,
                             show_support=getattr(args,'show_support',70.0),keep_tree_order=getattr(args,'keep_tree_order',False))
    print("\n[SUCCESS] ACVirus Tree-Pro completed!")

if __name__=="__main__":
    main()
