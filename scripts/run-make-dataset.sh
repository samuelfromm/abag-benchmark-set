#!/bin/bash

AADaM_script="/home/sfromm/git/AADaM-fork/AntibodyAntigenDatasetMaker.py"
input_db="/home/sfromm/git/abag-benchmark-set/all_structures/raw"
outdir="/home/sfromm/git/abag-benchmark-set/processed"

std_flags="
--inputDb=${input_db}
--date=2021/09/30
--abCompSeqCut=80
--withinDatasetCut=80
--outDir=${outdir}
--resCut=3.5
--cutoffStrict=True
--nx=True
--globalSeqID=True
--minAtomSeqresFraction=0.9
"

# parser.add_argument("-i", "--inputDb", dest="inputDb", help="path to the input database, which should be a complete download of all antibodies in complex with protein antigens, and have a summary file ending in 'summary.tsv'")
# parser.add_argument("-d", "--date", dest="date0", help="the date on or after which you want to make your test set; format YYYY/MM/DD")
# parser.add_argument("-cd", "--cutoffDate", dest="cutoffDate", help="a date past which you won't accept structures, useful for making training sets; format YYYY/MM/DD")
# parser.add_argument("-c1", "--abCompSeqCut", dest="abCompSeqCut", help="percent sequence ID minimum to disallow an ab-ag complex, based on past ab-ag complexes, comparing H loops to H loops and L loops to L loops")
# parser.add_argument("-c2", "--withinDatasetCut", dest="withinDatasetCut", help="percent sequence ID minimum used for knocking out complexes that are too similar w/ in the database, based on H or L or any antigen chain to any other antigen chain sequence ID")
# parser.add_argument("-m", "--methodsAllowed", dest="methodsAllowed", help="methods allowed, separated by a comma, like 'X-RAY DIFFRACTION,ELECTRON MICROSCOPY'")
# parser.add_argument("-w", "--whiteList", dest="whiteList", help="a comma separated list of pdb files that you want to definitely be included in the dataset (will still be filtered by resolution, sequence identity, method, etc. - this just ensures when structures 'knock each other out' the whitelisted ones will be given preference to remain. Useful for adding on new structures to a previously-made dataset)'")
# parser.add_argument("-o", "--outDir", dest="outDbPath", help="")
# parser.add_argument("-r", "--resCut", dest="resCut", help="defaults to 100")
# parser.add_argument("-cs", "--cutoffStrict", dest="cutoffStrict", help="make the cutoffs for sequence ID within the group strict, by applying them to the antigen and antibody loop sequences individually, instead of in combination (i.e. strict cutoffs would not allow an antibody-antigen complex with the same antigen as a previously accepted structure but different antibody in, while standard / not strict cutoffs would allow it)")
# parser.add_argument("-nx", "--nx", dest="nx", help="strip unnatural residues from the ends, and don't allow sequences with unnatural residues in the middle")
# parser.add_argument("-g", "--globalSeqID", dest="globalSeqID", help="use global sequence ID, number of matches divided by the length of the shorter sequence, to get sequence ID %. This can make the cutoff more stringent in practice, as a shorter sequence can map well to a much larger non-similar one by allowing many gaps. The default (when this flag is not given) is to do local alignment with gap penalites, then apply the sequence identity cutoffs to the best local alignment, if it's over 30 residues long.")

python3 ${AADaM_script} ${std_flags} --methodsAllowed='X-RAY DIFFRACTION,ELECTRON MICROSCOPY'