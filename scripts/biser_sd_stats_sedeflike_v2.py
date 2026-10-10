#!/usr/bin/env python3
"""Summarize a final SD callset from filter_biser_sd_sedeflike_v2.py.

For comparison with Vollger et al. (Science 2022), use
``merged_genomic_sd_span_bp`` / ``merged_genomic_sd_span_Mb``. This is the
union of the genomic coordinates of *both* SD copies, mirroring the paper's
SEDEF summary code (symmetric copies followed by BED merge). Do not use the
sum of pairwise alignment lengths as genome-wide SD content.
"""
from __future__ import annotations
import argparse, csv, gzip, re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

BINS=[(1000,2000,'1-2kb'),(2000,5000,'2-5kb'),(5000,10000,'5-10kb'),(10000,50000,'10-50kb'),(50000,100000,'50-100kb'),(100000,500000,'100-500kb'),(500000,None,'>=500kb')]
REQ={'chrom1','start1','end1','chrom2','start2','end2','sd_id','strand1','strand2','aligned_bp','aln_len'}

def nk(s): return [int(x) if x.isdigit() else x.lower() for x in re.split(r'(\d+)',s)]
def open_text(p): return gzip.open(p,'rt') if str(p).endswith('.gz') else open(p)
def merge(xs):
    xs=sorted(xs)
    out=[]
    for s,e in xs:
        if not out or s>out[-1][1]: out.append([s,e])
        else: out[-1][1]=max(out[-1][1],e)
    return [tuple(x) for x in out]
def bname(x):
    for lo,hi,n in BINS:
        if x>=lo and (hi is None or x<hi): return n
    return '<1kb'
def write_tsv(p,header,rows):
    with open(p,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(header); w.writerows(rows)
def read_fai(p):
    ans={}
    with open(p) as f:
        for l in f:
            if l.strip():
                a=l.rstrip().split('\t'); ans[a[0]]=int(a[1])
    return ans

def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-i','--input',required=True)
    ap.add_argument('-o','--out-prefix',required=True)
    ap.add_argument('--fai')
    args=ap.parse_args()
    calls=[]
    with open_text(args.input) as f:
        rd=csv.DictReader(f,delimiter='\t')
        if rd.fieldnames is None or not REQ <= set(rd.fieldnames):
            missing=', '.join(sorted(REQ-set(rd.fieldnames or [])))
            raise SystemExit('ERROR: not a sedeflike-v2 final table; missing '+missing)
        for n,r in enumerate(rd,2):
            for k in ('start1','end1','start2','end2','aligned_bp','aln_len'):
                r[k]=int(float(r[k]))
            if r['end1']<=r['start1'] or r['end2']<=r['start2']:
                raise SystemExit(f'ERROR row {n}: non-positive interval')
            r['max_mate_span_bp']=max(r['end1']-r['start1'],r['end2']-r['start2'])
            r['orientation']='direct' if r['strand1']==r['strand2'] else 'inverted'
            r['location']='intra' if r['chrom1']==r['chrom2'] else 'inter'
            calls.append(r)
    if not calls: raise SystemExit('ERROR: no calls')
    calls.sort(key=lambda r:(nk(r['chrom1']),r['start1'],nk(r['chrom2']),r['start2']))
    loci=defaultdict(list); copies={}
    orientation=Counter(); pair_aln=Counter(); pair_maxspan=Counter()
    for r in calls:
        orientation[(r['location'],r['orientation'])]+=1
        pair_aln[(r['location'],r['orientation'])]+=r['aln_len']
        pair_maxspan[(r['location'],r['orientation'])]+=r['max_mate_span_bp']
        for side in (1,2):
            c,s,e,st=r[f'chrom{side}'],r[f'start{side}'],r[f'end{side}'],r[f'strand{side}']
            loci[c].append((s,e)); copies[(c,s,e,st)]=1
    merged={c:merge(xs) for c,xs in loci.items()}
    total=sum(e-s for xs in merged.values() for s,e in xs)
    prefix=args.out_prefix
    with open(prefix+'.mates.bed','w') as out:
        for i,(c,s,e,st) in enumerate(sorted(copies,key=lambda x:(nk(x[0]),x[1],x[2],x[3])),1):
            out.write(f'{c}\t{s}\t{e}\tSDcopy{i:07d}\t0\t{st}\n')
    with open(prefix+'.merged_loci.bed','w') as out:
        i=1
        for c in sorted(merged,key=nk):
            for s,e in merged[c]:
                out.write(f'{c}\t{s}\t{e}\tSDlocus{i:07d}\n'); i+=1
    for metric,vals,fn in [('aln_len', [r['aln_len'] for r in calls], '.length_distribution.aln_len.tsv'),('max_mate_coordinate_span', [r['max_mate_span_bp'] for r in calls], '.length_distribution.max_mate_span.tsv')]:
        cnt=Counter(bname(x) for x in vals); bp=Counter()
        for x in vals: bp[bname(x)]+=x
        rows=[]
        for lo,hi,lab in BINS:
            rows.append([lab,lo,'' if hi is None else hi,cnt[lab],bp[lab],round(cnt[lab]/len(vals)*100,6),round(bp[lab]/sum(vals)*100,6)])
        write_tsv(prefix+fn,['length_definition','length_bin','lower_bp_inclusive','upper_bp_exclusive','n_sd_pairs','sum_bp','pair_fraction_pct','bp_fraction_pct'],[[metric]+r for r in rows])
    orows=[]
    for loc in ('intra','inter'):
        for ori in ('direct','inverted'):
            n=orientation[(loc,ori)]
            orows.append([loc,ori,n,round(100*n/len(calls),6),pair_aln[(loc,ori)],pair_maxspan[(loc,ori)]])
    write_tsv(prefix+'.orientation_summary.tsv',['chromosomal_class','orientation_class','n_sd_pairs','pair_fraction_pct','sum_aln_len_bp','sum_max_mate_span_bp'],orows)
    fai=read_fai(args.fai) if args.fai else {}
    chroms=set(merged)|set(fai)
    crows=[]
    for c in sorted(chroms,key=nk):
        span=sum(e-s for s,e in merged.get(c,[])); clen=fai.get(c,'')
        cov=round(100*span/clen,6) if isinstance(clen,int) and clen else ''
        ncopy=sum(1 for x in copies if x[0]==c)
        crows.append([c,clen,ncopy,len(merged.get(c,[])),span,cov])
    write_tsv(prefix+'.chromosome_distribution.tsv',['chrom','chromosome_length_bp','n_exact_sd_copy_intervals','n_merged_sd_loci','merged_sd_span_bp','merged_sd_coverage_pct'],crows)
    rows=[
        ['count','final_sd_pair_calls',len(calls)],
        ['count','exact_sd_copy_intervals',len(copies)],
        ['count','merged_sd_loci',sum(len(x) for x in merged.values())],
        ['length','sum_pair_alignment_span_bp',sum(r['aln_len'] for r in calls)],
        ['length','median_pair_alignment_span_bp',round(median(r['aln_len'] for r in calls),2)],
        ['span','sum_of_all_copy_coordinate_spans_bp',sum((r['end1']-r['start1'])+(r['end2']-r['start2']) for r in calls)],
        ['span','merged_genomic_sd_span_bp',total],
        ['span','merged_genomic_sd_span_Mb',round(total/1e6,6)],
    ]
    if fai:
        gl=sum(fai.values()); rows += [['reference','genome_length_bp',gl],['reference','merged_sd_coverage_pct',round(100*total/gl,6)]]
    write_tsv(prefix+'.summary.tsv',['section','metric','value'],rows)
    print(f'Done. {len(calls):,} pairs; merged genomic SD span = {total/1e6:.3f} Mb')

if __name__=='__main__': main()
