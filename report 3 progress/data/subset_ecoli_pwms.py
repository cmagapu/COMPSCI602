# subset_ecoli_pwms.py
ecoli_tfs = {
    "CRP","FNR","LexA","PurR","Lrp","IHF","Fis",
    "NarL","NarP","ArcA","PhoP","OmpR","OxyR",
    "SoxS","LacI","Ada","Cra","CpxR","Hns","LeuO"
}

infile  = "tf_pwms.meme"
outfile = "ecoli_pwms.meme"

with open(infile) as fin, open(outfile, "w") as fout:
    keep = False
    for line in fin:
        if line.startswith("MOTIF "):
            tfname = line.strip().split()[-1]
            keep = (tfname in ecoli_tfs)
        if keep:
            fout.write(line)

