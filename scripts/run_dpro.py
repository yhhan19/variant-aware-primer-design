import subprocess
from random import randint, seed
import primer3
import re
import matplotlib.pyplot as plt
import numpy as np
import time
from olivar.tiling_helper import generate_context, design_context_seq, get_primer, optimize
from olivar import build
from olivar.design import PrimerSetBadnessFast
import os
import tracemalloc

def read_fasta(path):
    sequences = {}
    with open(path, "r") as f:
        header = None
        seq_chunks = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):           # new record
                if header:
                    sequences[header] = ("".join(seq_chunks)).upper()
                header = line[1:]              # drop ">"
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header:
            sequences[header] = ("".join(seq_chunks)).upper()
    return sequences

def write_fasta(file_name, seqs):
    with open(file_name, "w") as fout:
        for label, seq in seqs.items():
            fout.write(">" + label + "\n")
            for i in range(0, len(seq), 80):
                fout.write(seq[i : i + 80] + "\n")


def rc(s):
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(comp)[::-1]

def upper(s):
    return re.sub(r'[^ACGTacgt]', '', s).upper()

def optimize_PDR(len_PDR, risk, ref, opt = "risk-d1", pref = ""):
    print("optimizing PDRs ...")
    size = len(ref)
    risk_optimizer_input = pref + ".risk_optimizer_input.txt"
    risk_optimizer_output = pref + ".risk_optimizer_output.txt"
    with open(risk_optimizer_input, "w") as file:
        for i in range(size):
            file.write(str(i) + "," + ref[i] + "," + str(risk[i]) + "\n")
    result = subprocess.run(
        ["./risk.exe", 
            "-i", risk_optimizer_input, 
            "-o", risk_optimizer_output,
            "-x", opt,
            "-S", "12345",
            "-Lp", str(len_PDR)
        ])
    with open(risk_optimizer_output, "r") as file:
        line = file.readline()
        words = line.split()
        PDRs = [[int(words[i]), int(words[i]) + len_PDR, int(words[i+1]), int(words[i+1]) + len_PDR] for i in range(len(words)) if i % 2 == 0]
    print("done")
    return PDRs

def optimize_dimer(results, pref = ""):
    print("optimizing dimers ...")
    start = time.perf_counter()
    K = len(results)
    N = max([len(result) for result in results])
    hairpin = []
    heterodimer = []
    for i in range(K):
        for j in range(N):
            index = i * N + j
            if j >= len(results[i]) or results[i][j] == "":
                hairpin.append([index, -1e9])
            else:
                hairpin.append([index, primer3_hairpin(results[i][j])])
            for i_ in range(i + 1, K):
                for j_ in range(N):
                    index_ = i_ * N + j_
                    if j >= len(results[i]) or j_ >= len(results[i_]) or results[i][j] == "" or results[i_][j_] == "":
                        heterodimer.append([index, index_, -1e9])
                    else:
                        heterodimer.append([index, index_,  primer3_dimer(results[i][j], results[i_][j_])])
    dimer_optimizer_input = pref + ".dimer_optimizer_input.txt"
    dimer_optimizer_output = pref + ".dimer_optimizer_output.txt"
    with open(dimer_optimizer_input, "w") as fout:
        fout.write(str(K) + " " + str(N) + " " + str(len(heterodimer)) + "\n")
        for h in hairpin:
            fout.write(str(h[0]) + " " + str(h[1]) + "\n")
        for h in heterodimer:
            fout.write(str(h[0]) + " " + str(h[1]) + " " + str(h[2]) + "\n")
    result = subprocess.run(
        ["./risk.exe", 
            "-i", dimer_optimizer_input, 
            "-o", dimer_optimizer_output,
            "-x", "dimer-h",
            "-S", "12345"
            "-I", "10000"
        ])
    with open(dimer_optimizer_output, "r") as fin:
        line = fin.readline()
        indices = [int(word) for word in line.split()]
    for i in range(0, K, 2):
        print(int(i / 2) + 1, results[i][indices[i]], results[i+1][indices[i+1]])
    primers_indices = []
    primers = []
    for i in range(0, K):
        primers.append(results[i][indices[i]])
        primers_indices.append(indices[i])
    print("done")
    end = time.perf_counter()
    print(end - start, "seconds")
    return primers, primers_indices

def optimize_dimer_badness(results, self_badness, pref = ""):
    print("optimizing dimers ...")
    start = time.perf_counter()
    K = len(results)
    N = max([len(result) for result in results])
    hairpin = []
    heterodimer = []
    for i in range(K):
        for j in range(N):
            index = i * N + j
            if j >= len(results[i]) or results[i][j] == "":
                hairpin.append([index, -1e9])
            else:
                hairpin.append([index, - self_badness[i][j] - two_primer_badness(results[i][j], results[i][j])])
            for i_ in range(i + 1, K):
                for j_ in range(N):
                    index_ = i_ * N + j_
                    if j >= len(results[i]) or j_ >= len(results[i_]) or results[i][j] == "" or results[i_][j_] == "":
                        heterodimer.append([index, index_, -1e9])
                    else:
                        heterodimer.append([index, index_, - two_primer_badness(results[i][j], results[i_][j_]) - two_primer_badness(results[i_][j_], results[i][j])])
    dimer_optimizer_input = pref + ".dimer_optimizer_input.txt"
    dimer_optimizer_output = pref + ".dimer_optimizer_output.txt"
    with open(dimer_optimizer_input, "w") as fout:
        fout.write(str(K) + " " + str(N) + " " + str(len(heterodimer)) + "\n")
        for h in hairpin:
            fout.write(str(h[0]) + " " + str(h[1]) + "\n")
        for h in heterodimer:
            fout.write(str(h[0]) + " " + str(h[1]) + " " + str(h[2]) + "\n")
    result = subprocess.run(
        ["./risk.exe", 
            "-i", dimer_optimizer_input, 
            "-o", dimer_optimizer_output,
            "-x", "dimer-h",
            "-S", "12345",
            "-I", "10000"
        ])
    with open(dimer_optimizer_output, "r") as fin:
        line = fin.readline()
        indices = [int(word) for word in line.split()]
    for i in range(0, K, 2):
        print(int(i / 2) + 1, results[i][indices[i]], results[i+1][indices[i+1]])
    primers = []
    badness = []
    for i in range(0, K):
        primers.append(results[i][indices[i]])
        badness.append(self_badness[i][indices[i]])
    print("done")
    end = time.perf_counter()
    print(end - start, "seconds")
    return primers, badness

def primer3_dimer(a, b):
    dg = primer3.bindings.calc_heterodimer(
        a,
        b,
        mv_conc=180.0,   # 0.18 M from Olivar -> 180 mM
        temp_c=60.0,     # Olivar annealing temperature
    ).dg + primer3.bindings.calc_heterodimer(
        b,
        a,
        mv_conc=180.0,   # 0.18 M from Olivar -> 180 mM
        temp_c=60.0,     # Olivar annealing temperature
    ).dg
    return dg / 2

def primer3_hairpin(a):
    dg = primer3.calc_hairpin(
        a,
        mv_conc=180.0,   # 0.18 M from Olivar -> 180 mM
        temp_c=60.0,     # Olivar annealing temperature
    ).dg
    return dg
    
def evaluate(primers):
    final_hairpin = []
    final_heterodimer = []
    for i in range(0, len(primers)):
        final_hairpin.append(primer3_hairpin(primers[i]))
        for j in range(i + 1, len(primers)):
            final_heterodimer.append(primer3_dimer(primers[i], primers[j]))
    print(len(final_hairpin), np.mean(final_hairpin), np.std(final_hairpin))
    print(len(final_heterodimer), np.mean(final_heterodimer), np.std(final_heterodimer))
    return final_hairpin, final_heterodimer

def evaluate_badness(primers, self_badness):
    final_hairpin = []
    final_heterodimer = []
    for i in range(0, len(primers)):
        final_hairpin.append(self_badness[i] + two_primer_badness(primers[i], primers[i]))
        for j in range(i + 1, len(primers)):
            final_heterodimer.append(two_primer_badness(primers[i], primers[j]) + two_primer_badness(primers[j], primers[i]))
    print(len(final_hairpin), np.mean(final_hairpin), np.std(final_hairpin))
    print(len(final_heterodimer), np.mean(final_heterodimer), np.std(final_heterodimer))
    print(sum(final_hairpin) + sum(final_heterodimer))
    return final_hairpin, final_heterodimer + final_hairpin

def plot_histogram(axes, data, color, title, xlabel, x = None, y = None, invert=True):
    axes.set_title(title, fontsize=18)
    axes.set_xlabel(xlabel, fontsize=14)
    axes.set_ylabel('frequency', fontsize=14)
    if y != None:
        axes.set_ylim(y[0], y[1])
        mean = np.mean(data)
        std = np.std(data)
        axes.axvline(mean, color='red', linestyle='dashed', linewidth=2)
        axes.text(mean + 0.2, axes.get_ylim()[1]*0.9, f'Mean = {mean:.2f}', color='red')
        # axes.axvline(mean - std, color='orange', linestyle='-', linewidth=0.5, label=f'Std = {std:.2f}')
        # axes.axvline(mean + std, color='orange', linestyle='-', linewidth=0.5)
        # axes.text(mean - std - 0.2, axes.get_ylim()[1]*0.8, f'Std = {std:.2f}', color='orange')
    if invert: 
        axes.set_xlim(x[0], x[1])
        bin_width = (x[1] - x[0]) / 50
        bins = np.arange(x[0], x[1] + bin_width, bin_width)
        axes.hist(data, bins=bins, color=color, edgecolor='black')
        axes.invert_xaxis()
    else:
        bins = np.logspace(np.log10(1e-2), np.log10(100), 50)
        axes.hist(data, bins=bins, color=color, edgecolor='black')
        axes.set_xscale("log")

def plot_all(data, output):
    x = [
            [
                [min(data[0][0] + data[0][1]), max(data[0][0] + data[0][1])],
                [min(data[0][0] + data[0][1]), max(data[0][0] + data[0][1])],
                [min(data[0][2] + data[0][3]), max(data[0][2] + data[0][3])],
                [min(data[0][2] + data[0][3]), max(data[0][2] + data[0][3])],
            ],
            [
                [min(data[1][0] + data[1][1]), max(data[1][0] + data[1][1])],
                [min(data[1][0] + data[1][1]), max(data[1][0] + data[1][1])],
                [min(data[1][2] + data[1][3]), max(data[1][2] + data[1][3])],
                [min(data[1][2] + data[1][3]), max(data[1][2] + data[1][3])],
            ],
        ]
    print(x)
    y = [
            [
                [0, -1e9],
                [0, -1e9],
                [0, -1e9],
                [0, -1e9],
            ],
            [
                [0, -1e9],
                [0, -1e9],
                [0, -1e9],
                [0, -1e9],
            ],
        ]
    colors = [
        ['skyblue', 'salmon', 'skyblue', 'salmon'],
        ['lightgreen', 'plum', 'lightgreen', 'plum']
    ]
    titles = [
        ['Badness$_{SADDLE}$', 'Badness$_{\Delta-PRO}$', 'Badness$_{SADDLE}$', 'Badness$_{\Delta-PRO}$'],
        ['Heterodimer$_{SADDLE}$', 'Heterodimer$_{\Delta-PRO}$', 'Heterodimer$_{SADDLE}$', 'Heterodimer$_{\Delta-PRO}$'],
    ]
    xlabels = ['Badness value', '$\Delta G$ value']
    fig, axes = plt.subplots(2, 4, figsize=(16, 10))
    for i in range(2):
        for j in range(4):
            plot_histogram(axes[i][j], data[i][j], colors[i][j], titles[i][j], '', x[i][j], invert=i==1)
            y[i][j][1] = axes[i][j].get_ylim()[1]
            y[i][j][0] = axes[i][j].get_ylim()[0]
    y = [
            [
                [min(y[0][0][0], y[0][1][0]), max(y[0][0][1], y[0][1][1])],
                [min(y[0][0][0], y[0][1][0]), max(y[0][0][1], y[0][1][1])],
                [min(y[0][2][0], y[0][3][0]), max(y[0][2][1], y[0][3][1])],
                [min(y[0][2][0], y[0][3][0]), max(y[0][2][1], y[0][3][1])],
            ],
            [
                [min(y[1][0][0], y[1][1][0]), max(y[1][0][1], y[1][1][1])],
                [min(y[1][0][0], y[1][1][0]), max(y[1][0][1], y[1][1][1])],
                [min(y[1][2][0], y[1][3][0]), max(y[1][2][1], y[1][3][1])],
                [min(y[1][2][0], y[1][3][0]), max(y[1][2][1], y[1][3][1])],
            ],
        ]
    plt.close(fig)
    fig, axes = plt.subplots(2, 4, figsize=(20, 6),  constrained_layout=True)
    for i in range(2):
        for j in range(4):
            plot_histogram(axes[i][j], data[i][j], colors[i][j], titles[i][j], xlabels[i], x[i][j], y[i][j], invert=i==1)
    fig.text(0.25, 0.96, "Pool 1", ha="center", va="bottom", fontsize=18)
    fig.text(0.73, 0.96, "Pool 2", ha="center", va="bottom", fontsize=18)
    fig.savefig(output + ".png")
    fig.savefig(output + ".pdf")

def olivar_optimize_PDR(risk, seed = 10):
    np.random.seed(seed)
    size = len(risk)
    N = 500 * size // 420
    rand_int = np.random.randint(0, 2**32, size = N)
    design = [generate_context((risk, 1, size, 420, 252, rand_int[i])) for i in range(N)]
    all_loss = [d[2] for d in design]
    all_loss.sort(reverse=True)
    best_design = sorted(design, key=lambda x:x[2])
    all_context_seq, all_risk, loss = best_design[0]
    print('Loss of the best design: %.3f' % loss)
    return best_design[0][0].tolist()

def get_concensus(seqs):
    concensus = ""
    lens = [len(seq) for _, seq in seqs.items()]
    for i in range(lens[0]):
        count = {}
        count['-'] = 0
        for _, seq in seqs.items():
            if not seq[i] in count:
                count[seq[i]] = 0
            count[seq[i]] += 1
        c = max([[key, value] for key, value in count.items() if key in ['A', 'C', 'G', 'T']], key=lambda kv: kv[1])
        concensus += c[0]
    return concensus

def tiling_copied(ref_path: str, out_path: str, title: str, max_amp_len: int, min_amp_len: int, 
    w_egc: float, w_lc: float, w_ns: float, w_var: float, w_sensi: float, w_combi: float,
    temperature: float, salinity: float, dG_max: float, min_GC: float, max_GC: float, 
    min_complexity: float, max_len: int, check_var: bool, fP_prefix: str, rP_prefix: str,
    seed: int, threads: int, iterMul: int, deg: bool):
    '''
    Design tiled amplicons. 
    Input:
        ref_path: Path to the Olivar reference file (.olvr), or the directory of reference files for multiple targets.
        out_path: Output directory [./].
        title: Name of this design [olivar-design].
        max_amp_len: Maximum amplicon length [420].
        min_amp_len: Minimum amplicon length. 0.9*{max-amp-len} if set as None. Minimum 120.
        w_egc: Weight for extreme GC content [1.0].
        w_lc: Weight for low sequence complexity [1.0].
        w_ns: Weight for non-specificity [1.0].
        w_var: Weight for variations [1.0].
        
        w_sensi: Weight for sensitivity [1.0].
        w_combi: Weight for combinations [1.0].
        
        temperature: PCR annealing temperature [60.0].
        salinity: Concentration of monovalent ions in units of molar [0.18].
        dG_max: Maximum free energy change of a primer in kcal/mol [-11.8].
        min_GC: Minimum GC content of a primer [0.2].
        max_GC: Maximum GC content of a primer [0.75].
        min_complexity: Minimum sequence complexity of a primer [0.4].
        max_len: Maximum length of a primer [36].
        check_var: Filter out primer candidates with variations within 5nt of 3' end [False]. 
            Setting check_var=True is not recommended when a lot of variations are provided, 
            since this would significantly reduce the number of primer candidates. 
        fP_prefix: Prefix of forward primer. Empty string '' by default.
        rP_prefix: Prefix of reverse primer. Empty string '' by default.
        seed: Random seed for optimizing primer design regions and primer dimer [10].
        threads: Number of threads [1].
        iterMul: Multiplier of iterations during PDR optimization [1].
        deg: Control whether use degenerate mode or not.
    '''
    
    config = {
        'ref_path': ref_path, 
        'out_path': out_path, 
        'title': title, 
        'max_amp_len': max_amp_len, 
        'min_amp_len': min_amp_len, 
        'w_egc': w_egc, 
        'w_lc': w_lc, 
        'w_ns': w_ns, 
        'w_var': w_var, 
        
        'w_sensi': w_sensi,
        'w_combi': w_combi,
        
        'temperature': temperature, 
        'salinity': salinity, 
        'dG_max': dG_max, 
        'min_GC': min_GC, 
        'max_GC': max_GC, 
        'min_complexity': min_complexity, 
        'max_len': max_len, 
        'check_SNP': check_var, 
        'fP_prefix': fP_prefix, 
        'rP_prefix': rP_prefix, 
        'seed': seed, 
        'threads': threads,
        'iterMul': iterMul
    }
    print("copy")
    if not os.path.exists(config['out_path']):
        os.makedirs(config['out_path'])

    # store config values
    design_ref_path = config['ref_path']

    # validate input .olvr file(s)
    ref_path_dict = dict() # {ref_name: ref_path}

    REFEXT = '.olvr'
    # check input is a file or a directory
    if os.path.isfile(design_ref_path) and design_ref_path.endswith(REFEXT):
        # use the name of the .olvr file as reference name
        file = os.path.basename(design_ref_path)
        ref_path_dict[file[:-len(REFEXT)]] = design_ref_path
    elif os.path.isdir(design_ref_path):
        for file in sorted(os.listdir(design_ref_path)):
            file_path = os.path.join(design_ref_path, file)
            if os.path.isfile(file_path) and file_path.endswith(REFEXT):
                # use file name as reference name
                ref_path_dict[file[:-len(REFEXT)]] = file_path
        if ref_path_dict:
            1
        else:
            raise FileNotFoundError(f'No {REFEXT} file found in the directory "{design_ref_path}".')
    else:
        raise FileNotFoundError(f'Input is neither a {REFEXT} file nor a directory.')
    
    # design PDRs for each reference
    all_plex_info = {}
    all_ref_info = {}
    for ref_name, ref_path in ref_path_dict.items():
        config['ref_path'] = ref_path
        temp_dict, risk_arr, gc_arr, comp_arr, hits_arr, var_arr, sensi_arr, combi_arr, all_loss, seq_record = design_context_seq(config, deg=deg)
        all_plex_info.update(temp_dict)
        all_ref_info[ref_name] = {
            'risk_arr': risk_arr, 
            'gc_arr': gc_arr, 
            'comp_arr': comp_arr, 
            'hits_arr': hits_arr, 
            'var_arr': var_arr, 

            'sensi_arr': sensi_arr,
            'combi_arr': combi_arr,

            'all_loss': all_loss, 
            'seq_record': seq_record, 
        }
    
    # revert modified config values
    config['ref_path'] = design_ref_path
    all_plex_info_primer = get_primer(all_plex_info, config)
    results = [[], []]
    badness = [[], []]
    for plex_id, plex_info in all_plex_info.items():
        fps = []
        fbs = []
        for i in range(len(plex_info['fP_candidate'])):
            fp = plex_info['fP_candidate'].iloc[i]['seq'][0]
            fps.append(fp)
            fbs.append(plex_info['fP_candidate'].iloc[i]['badness'])
        rps = []
        rbs = []
        for i in range(len(plex_info['rP_candidate'])):
            rp = plex_info['rP_candidate'].iloc[i]['seq'][0]
            rps.append(rp)
            rbs.append(plex_info['rP_candidate'].iloc[i]['badness'])
        results[plex_info['tube'] - 1].append(fps)
        results[plex_info['tube'] - 1].append(rps)
        badness[plex_info['tube'] - 1].append(fbs)
        badness[plex_info['tube'] - 1].append(rbs)

    all_plex_info_optimize, learning_curve = optimize(all_plex_info_primer, config)
    olivar_result = [[], []]
    olivar_badness = [[], []]
    for plex_id, plex_info in all_plex_info_optimize.items():
        fp = all_plex_info[plex_id]['fP']['seq'][0]
        rp = all_plex_info[plex_id]['rP']['seq'][0]
        olivar_result[plex_info['tube'] - 1].append(fp)
        olivar_result[plex_info['tube'] - 1].append(rp)
        olivar_badness[plex_info['tube'] - 1].append(all_plex_info[plex_id]['fP']['badness'])
        olivar_badness[plex_info['tube'] - 1].append(all_plex_info[plex_id]['rP']['badness'])
        print(plex_info['tube'], fp, rp)
    return results, olivar_result, badness, olivar_badness

def get_rates_substitute(seqs, ref):
    rates = []
    lens = [len(seq) for _, seq in seqs.items()]
    for i in range(lens[0]):
        count = {}
        count['-'] = 0
        for _, seq in seqs.items():
            if not seq[i] in count:
                count[seq[i]] = 0
            count[seq[i]] += 1
        total = sum([value for key, value in count.items() if key in ['A', 'C', 'G', 'T']])
        rates.append(sum([value for key, value in count.items() if key != ref[i] and key in ['A', 'C', 'G', 'T']]) / total + count['-'] / len(seqs.items()))
    return rates

def evaluate_olivar(primers, self_badness):
    fp = [primers[i] for i in range(len(primers)) if i % 2 ==0]
    rp = [primers[i] for i in range(len(primers)) if i % 2 ==1]
    return PrimerSetBadnessFast(
        all_fP=[fp], 
        all_rP=[rp],
        existing=[],
    )[0] + sum([b for b in self_badness])

from olivar.basic import revcomp
from collections import defaultdict

def pair_badness(p, b):
    sb = sum(b)
    sp, _ = PrimerSetBadnessFast(all_fP = p, all_rP = [], existing=[])
    print(sb, sp, sp + sb)
    total = 0
    for i in range(len(p)):
        for j in range(len(p)):
            total += two_primer_badness(p[i], p[j])
    print(total, total + sb)

GC_LETTER = set(["C", "G", "c", "g"])
PENALTY_OFFSET = 0
END4 = 1
END5 = 4
END6 = 20
MIDDLE7 = 100
MIDDLE8 = 500

def hash_weights_for_primer(primer: str, conc: float = 1.0):
    """
    Compute per-primer hash weights w_hash(type,kmer) exactly like SADDLE's
    endhash4/5/6 and middlehash7/8 building step.

    Returns a dict: {("end4", kmer): weight, ..., ("mid8", kmer): weight}
    """
    p = primer.lower()
    l = len(p)
    w = defaultdict(float)

    # End 4/5/6 (last bases only)
    if l >= 4:
        kmer = p[-4:]
        w[("end4", kmer)] += conc
    if l >= 5:
        kmer = p[-5:]
        w[("end5", kmer)] += conc
    if l >= 6:
        kmer = p[-6:]
        w[("end6", kmer)] += conc

    # Middle 7-mers
    if l >= 7:
        for j in range(l - 6):
            kmer = p[j:j+7]
            # distance from end of kmer to 3' end: end index = j+6, 3' index = l-1
            d_other = (l - 1) - (j + 6)  # = l - j - 7
            weight = conc * (PENALTY_OFFSET + 1) / (d_other + 1 + PENALTY_OFFSET)
            w[("mid7", kmer)] += weight

    # Middle 8-mers
    if l >= 8:
        for j in range(l - 7):
            kmer = p[j:j+8]
            d_other = (l - 1) - (j + 7)  # = l - j - 8
            weight = conc * (PENALTY_OFFSET + 1) / (d_other + 1 + PENALTY_OFFSET)
            w[("mid8", kmer)] += weight

    return w

def query_weights_for_primer(primer: str, conc: float = 1.0):
    """
    Compute per-primer query weights w_query(type,kmer) from the reverse
    complement, matching the (PENALTY_OFFSET+1)/(j+1+PENALTY_OFFSET) factor
    used when scanning one_side_primer in SADDLE.

    Returns a dict: {("end4", kmer): sum_j w_query, ..., ("mid8", kmer): sum_j w_query}
    """
    p = primer.lower()
    c = revcomp(p)
    l = len(c)
    wq = defaultdict(float)

    # 4-mers (end-like, with END4 coefficient later)
    for j in range(l - 3):
        kmer = c[j:j+4]
        w_query = conc * (PENALTY_OFFSET + 1) / (j + 1 + PENALTY_OFFSET)
        wq[("end4", kmer)] += w_query

    # 5-mers
    if l >= 5:
        for j in range(l - 4):
            kmer = c[j:j+5]
            w_query = conc * (PENALTY_OFFSET + 1) / (j + 1 + PENALTY_OFFSET)
            wq[("end5", kmer)] += w_query

    # 6-mers
    if l >= 6:
        for j in range(l - 5):
            kmer = c[j:j+6]
            w_query = conc * (PENALTY_OFFSET + 1) / (j + 1 + PENALTY_OFFSET)
            wq[("end6", kmer)] += w_query

    # 7-mers (middle)
    if l >= 7:
        for j in range(l - 6):
            kmer = c[j:j+7]
            w_query = conc * (PENALTY_OFFSET + 1) / (j + 1 + PENALTY_OFFSET)
            wq[("mid7", kmer)] += w_query

    # 8-mers (middle)
    if l >= 8:
        for j in range(l - 7):
            kmer = c[j:j+8]
            w_query = conc * (PENALTY_OFFSET + 1) / (j + 1 + PENALTY_OFFSET)
            wq[("mid8", kmer)] += w_query

    return wq

def two_primer_badness(a: str, b: str, conc_a: float = 1.0, conc_b: float = 1.0) -> float:
    """
    Badness(a, b) defined so that for any set S of primers,

        sum_{p in S} sum_{q in S} two_primer_badness(p, q)

    equals the SADDLE hash-based Badness L(S) (up to float noise),
    assuming same constants and concentration usage.
    """
    # Per-primer contributions (intrinsic to each primer)
    w_hash_b = hash_weights_for_primer(b, conc_b)      # "other" side (hash-like)
    w_query_a = query_weights_for_primer(a, conc_a)    # "query" side (reverse complement)

    total = 0.0

    for (kind, kmer), wq in w_query_a.items():
        wb = w_hash_b.get((kind, kmer), 0.0)
        if wb == 0.0:
            continue

        num_gc = sum(1 for ch in kmer if ch in GC_LETTER)

        if kind == "end4":
            coeff = END4
        elif kind == "end5":
            coeff = END5
        elif kind == "end6":
            coeff = END6
        elif kind == "mid7":
            coeff = MIDDLE7
        elif kind == "mid8":
            coeff = MIDDLE8
        else:
            continue  # should not happen

        # GC bonus 2^{numGC}, like PrimerSetBadnessFast
        factor = coeff * (2 ** num_gc)

        total += wq * wb * factor

    return total

def get_rates_gini(seqs):
    rates = []
    lens = [len(seq) for _, seq in seqs.items()]
    for i in range(lens[0]):
        count = {}
        count['-'] = 0
        for _, seq in seqs.items():
            if not seq[i] in count:
                count[seq[i]] = 0
            count[seq[i]] += 1
        total = sum([value for key, value in count.items() if key in ['A', 'C', 'G', 'T']])
        gini = 1 - sum([(value / total) ** 2 for key, value in count.items() if key in ['A', 'C', 'G', 'T']])
        rates.append(gini)
    return rates

def design(PDR, template, candidates):
    seq_args = {
        "SEQUENCE_ID": "example",
        "SEQUENCE_TEMPLATE": template,
        "SEQUENCE_INCLUDED_REGION": [PDR[0], PDR[3] - PDR[0]],
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [
            PDR[0], PDR[1] - PDR[0],
            PDR[2], PDR[3] - PDR[2],
        ]
    }
    global_args = {
        "PRIMER_NUM_RETURN": candidates,

        # Geometry — make sure your windows + this range are compatible
        "PRIMER_PRODUCT_SIZE_RANGE": [[60, 1200]],

        # Length & Tm — wide first, tighten later
        "PRIMER_MIN_SIZE": 16, "PRIMER_OPT_SIZE": 22, "PRIMER_MAX_SIZE": 34,
        "PRIMER_MIN_TM": 48.0, "PRIMER_OPT_TM": 60.0, "PRIMER_MAX_TM": 72.0,

        # Composition
        "PRIMER_MIN_GC": 30.0, "PRIMER_MAX_GC": 75.0,
        "PRIMER_MAX_POLY_X": 8,
        "PRIMER_MAX_NS_ACCEPTED": 2,

        # Thermodynamics / structure — loosen a lot
        "PRIMER_MAX_SELF_ANY_TH": 120.0,
        "PRIMER_MAX_SELF_END_TH": 120.0,
        "PRIMER_MAX_HAIRPIN_TH": 120.0,
        "PRIMER_MAX_END_STABILITY": 18.0,

        # Specificity (disable/ignore mispriming checks)
        # If you previously set a mispriming library, clear it or make thresholds huge:
        "PRIMER_MAX_TEMPLATE_MISPRIMING_TH": 1e9,
        "PRIMER_MAX_LIBRARY_MISPRIMING": 1e9,

        # Clamp off
        "PRIMER_GC_CLAMP": 0,

        # Thermo environment — lower DNA conc reduces Tm; raise if Tm is too low
        "PRIMER_DNA_CONC": 50.0,          # nM
        "PRIMER_SALT_MONOVALENT": 50.0,   # mM
        "PRIMER_SALT_DIVALENT": 1.5,      # mM
        "PRIMER_DNTP_CONC": 0.2,          # mM

        # Debug
        "PRIMER_EXPLAIN_FLAG": 1,
    }
    result = primer3.bindings.design_primers(seq_args, global_args)

    results = [[], []]
    results_pos = [[], []]
    for i in range(candidates):
        if not "PRIMER_LEFT_" + str(i) + "_SEQUENCE" in result:
            assert(i > 0)
            results[0].append("")
            results[1].append("")
            results_pos[0].append([-1, -1])
            results_pos[1].append([-1, -1])
            break
        left_seq  = result["PRIMER_LEFT_" + str(i) + "_SEQUENCE"]
        right_seq = result["PRIMER_RIGHT_" + str(i) + "_SEQUENCE"]
        size      = result["PRIMER_PAIR_" + str(i) + "_PRODUCT_SIZE"]
        left_pos  = result["PRIMER_LEFT_" + str(i)]
        right_pos = result["PRIMER_RIGHT_" + str(i)]
        amp_start = left_pos[0]
        amp_end   = right_pos[0] + right_pos[1]
        amplicon  = template[amp_start:amp_end]
        start = right_pos[0] - right_pos[1] + 1
        end_excl = right_pos[0] + 1
        footprint_forward = template[start:end_excl]
        print(f"Left  @ {left_pos} -> {left_seq}", template[left_pos[0]:left_pos[0]+left_pos[1]])
        print(f"Right @ {right_pos} -> {right_seq}", rc(footprint_forward), footprint_forward)
        print(f"Amplicon {amp_end-amp_start} bp: {amplicon}")
        het = primer3.calc_heterodimer(left_seq, right_seq)
        hl  = primer3.calc_hairpin(left_seq)
        hr  = primer3.calc_hairpin(right_seq)
        results[0].append(left_seq)
        results[1].append(right_seq)
        results_pos[0].append([left_pos[0], left_pos[0] + left_pos[1]])
        results_pos[1].append([right_pos[0] - right_pos[1] + 1, right_pos[0] + 1])
        print(template[left_pos[0]:left_pos[0] + left_pos[1]])
        print(template[right_pos[0] - right_pos[1] + 1:right_pos[0] + 1])
        assert(left_pos[0] >= PDR[0])
        assert(left_pos[0] + left_pos[1] <= PDR[1])
        assert(right_pos[0] - right_pos[1] + 1 >= PDR[2])
        assert(right_pos[0] + 1 <= PDR[3])

    region_id = f'size_{PDR[1] - PDR[0]}_region_{PDR[0]}_{PDR[2]}'
    primer_pairs = []
    num_pairs = result.get('PRIMER_PAIR_NUM_RETURNED', 0)
    for i in range(num_pairs):
        primer_info = {
            'pair_id': f"{region_id}_pair_{i}",
            'left_primer': {
                'sequence': result[f'PRIMER_LEFT_{i}_SEQUENCE'],
                'tm': result[f'PRIMER_LEFT_{i}_TM'],
                'gc_percent': result[f'PRIMER_LEFT_{i}_GC_PERCENT'],
                'position': result[f'PRIMER_LEFT_{i}'][0],
                'length': result[f'PRIMER_LEFT_{i}'][1]
            },
            'right_primer': {
                'sequence': result[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                'tm': result[f'PRIMER_RIGHT_{i}_TM'],
                'gc_percent': result[f'PRIMER_RIGHT_{i}_GC_PERCENT'],
                'position': result[f'PRIMER_RIGHT_{i}'][0],
                'length': result[f'PRIMER_RIGHT_{i}'][1]
            },
            'product_size': result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']
        }
        primer_pairs.append(primer_info)
    return results, region_id, primer_pairs, results_pos

def design_primers(candidates, ref):
    print("designing primers ...")
    results = [[], []]
    prism_results = [{}, {}]
    results_pos = [[], []]
    for i in range(len(PDRs)):
        print(i, PDRs[i])
        result, region_id, prism_result, result_pos = design(PDRs[i], ref, candidates)
        results[i % 2] += result
        prism_results[i % 2][region_id] = prism_result
        results_pos[i % 2] += result_pos
    print("done")
    return results, prism_results, results_pos

if __name__ == "__main__":

    seqs = read_fasta("./ZIKA.msa.fasta.filtered.fasta")
    concensus = get_concensus(seqs)
    rates = get_rates_gini(seqs)
    PDRs = optimize_PDR(40, rates, concensus)
    candidates = 10
    results, prism_results, results_pos = design_primers(candidates, concensus)
    primers = [optimize_dimer(results[i]) for i in [0, 1]]

    output = []
    for pool in [0, 1]:
        p = primers[pool]
        r = results[pool]
        ps = results_pos[pool]
        for i in range(len(p[0])):
            primer = p[0][i]
            index = p[1][i]
            primer_ = r[i][index]
            pos = ps[i][index]
            output.append([pool, i, pos[0], pos[1], primer])
    output = sorted(output, key=lambda x: x[2])
    with open("dpro.ZIKA.bed", "w") as f:
        for o in output:
            line = '\t'.join(['ref', str(o[2]), str(o[3]), 'dpro_' + str(o[1] // 2 * 2 + o[0]) + '_' + ['LEFT', 'RIGHT'][o[1] % 2], str(o[0] + 1), ['+', '-'][o[1] % 2], o[-1]])
            print(line)
            f.write(line + '\n')
