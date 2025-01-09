import numpy as np
import scipy.special
import scipy.stats
import pandas as pd
import tqdm
import torch
import fastDTWF
import os


# Generate wright_fisher_trajectories.txt
print('Generating wright_fisher_trajectories.txt')


def binom_stabilizing(N, k, s):
    p = k/N
    p00 = (1-p)**2
    p01 = 2*p*(1-p) * (1-s)
    p11 = p**2
    normalizer = (p00 + p01 + p11)
    p00 /= normalizer
    p01 /= normalizer
    p11 /= normalizer

    p_new = 0.5*p01 + p11
    return np.random.binomial(N, p_new)


np.random.seed(43)
trajectories = np.zeros((10000, 1000))
for i in tqdm.tqdm(range(10000)):
    while True:
        this_traj = [1]
        for k in range(999):
            this_traj.append(binom_stabilizing(20000, this_traj[-1], 0.001))
            # if this_traj[-1] == 0:
            #     break
        break
        # if len(this_traj) >= 1000 and this_traj[-1] != 0:
        #     break
    trajectories[i] = this_traj

np.savetxt('data/wright_fisher_trajectories.txt', trajectories)


# Generate stabilizing selection results
print('Running fastDTWF')

if (
    not os.path.isfile('data/s_hets_N_20000.npy')
    or not os.path.isfile('data/vgs_N_20000.npy')
    or not os.path.isfile('data/full_afs_N_20000.npy')
):
    N = 20000
    mu = torch.tensor(1.25e-8).double()
    s_hets = []
    vgs = []
    matrix = np.zeros((50, 2*N+1))
    for s_idx, s_het_np in tqdm.tqdm(
        enumerate(np.geomspace(1e-7, 0.05, num=50)), total=50
    ):
        s_hets.append(s_het_np)
        s_het = torch.tensor(s_het_np).double()
        s_hom = torch.tensor(0.).double()

        afs = fastDTWF.get_likelihood(
            pop_size_list=[2*N],
            switch_points=[0],
            sample_size=2*N,
            s_het=s_het,
            s_hom=s_hom,
            mu_0_to_1=mu,
            mu_1_to_0=mu*0,
            dtwf_tv_sd=0.1,
            dtwf_row_eps=1e-8,
            sampling_tv_sd=0.05,
            sampling_row_eps=1e-8,
            no_fix=True,
            sfs=False,
            injection_rate=0.,
        ).detach().numpy()

        matrix[s_idx] = afs

        freqs = np.arange(2*N+1) / (2*N)
        vgs.append(
            (2 * afs * freqs * (1-freqs)).sum()
        )

    np.save('data/s_hets_N_{}.npy'.format(N), s_hets)
    np.save('data/vgs_N_{}.npy'.format(N), vgs)
    np.save('data/full_afs_N_{}.npy'.format(N), matrix)


# Generate simulated_realized_heritability.txt
print('Generating simulated_realized_heritability.txt')
ss = np.load('data/s_hets_N_20000.npy')
vgs = np.load('data/vgs_N_20000.npy')
afs = np.load('data/full_afs_N_20000.npy')

np.random.seed(44)

x_samples = []
y_samples = []

sample_size = 50
max_f = 0
for s in tqdm.tqdm(np.linspace(ss[0] + 1e-10, ss[29]-1e-10, num=1000)):
    x_samples.extend([s]*sample_size)
    s_idx = np.searchsorted(ss, s)
    ps = (
        afs[s_idx] * (ss[s_idx+1] - s)
        / (ss[s_idx+1] - ss[s_idx])
        + afs[s_idx + 1] * (s - ss[s_idx])
        / (ss[s_idx+1] - ss[s_idx])
    )
    f = np.random.choice(40001, size=sample_size, replace=True, p=ps) / 40000
    f = np.random.binomial(600000, f) / 600000
    max_maf = np.max(np.minimum(f, 1-f))
    max_f = np.max([max_maf, max_f])
    y_samples.extend(2*f*(1-f)*ss[s_idx])

np.savetxt(
    'data/simulated_realized_heritability.txt',
    np.array([x_samples, y_samples]).T
)


# Generate p_val_ranking_sims.csv
print('Generating p_val_ranking_sims.csv')


def run_sims(N, p, sig_p, scale_factors):
    np.random.seed(42)
    ss = np.load('data/s_hets_N_20000.npy')
    afs = np.load('data/full_afs_N_20000.npy')
    num_traits = 18
    num_replicates = 100
    q1 = []
    q2 = []
    q3 = []
    q4 = []
    props_q1 = np.zeros((num_replicates, 4))
    props_q2 = np.zeros((num_replicates, 4))
    props_q3 = np.zeros((num_replicates, 4))
    props_q4 = np.zeros((num_replicates, 4))

    specificities = np.zeros((num_replicates, 4))
    frequencies = np.zeros((num_replicates, 4))
    scale_factors = [scale_factors] * num_replicates
    for i in tqdm.tqdm(range(num_replicates)):
        sf = scale_factors[i]
        alphas_sq = (1/sf)*np.exp(
            sf*scipy.stats.multivariate_normal(
                cov=(
                    p*np.eye(num_traits)*9
                    + (1-p)*np.ones((num_traits, num_traits))*9
                )
            ).rvs(10000000).T
        )*1e-7
        selection = alphas_sq.sum(axis=0)
        specificity = alphas_sq / selection
        sbin = np.maximum(np.searchsorted(ss, selection, 'right'), 0)
        freqs = np.zeros_like(sbin)
        for s in np.unique(sbin):
            count = (sbin == s).sum()
            if s == len(afs):
                continue
            freqs[sbin == s] = np.random.choice(40001, p=afs[s], size=count)

        alphas_sq = alphas_sq[:, freqs > 0]
        specificity = specificity[:, freqs > 0]
        freqs = freqs[freqs > 0]
        mafs = np.minimum(freqs / 40000, 1. - freqs/40000)
        test_stats = scipy.stats.ncx2(
            1, alphas_sq * 2 * freqs/40000 * (1-freqs/40000) * N
        ).rvs()
        pvals = scipy.stats.chi2(1).sf(test_stats)
        this_q1 = []
        this_q2 = []
        this_q3 = []
        this_q4 = []
        for j in range(num_traits):
            b1, b2, b3, b4 = np.percentile(
                pvals[j, pvals[j] < sig_p], [25, 50, 75, 100]
            )

            q1_keep = (pvals[j] >= 0)*(pvals[j] < b1)
            q1_num_traits = (pvals[:, q1_keep] < sig_p).sum(axis=0)
            this_q1.append((q1_num_traits-1).mean())
            props_q1[i, 0] += (q1_num_traits == 1).sum()
            props_q1[i, 1] += (q1_num_traits == 2).sum()
            props_q1[i, 2] += (q1_num_traits == 3).sum()
            props_q1[i, 3] += (q1_num_traits > 3).sum()

            specificities[i, 0] += specificity[j, q1_keep].mean() / num_traits
            frequencies[i, 0] += mafs[q1_keep].mean() / num_traits

            q2_keep = (pvals[j] >= b1)*(pvals[j] < b2)
            q2_num_traits = (pvals[:, q2_keep] < sig_p).sum(axis=0)
            this_q2.append((q2_num_traits-1).mean())
            props_q2[i, 0] += (q2_num_traits == 1).sum()
            props_q2[i, 1] += (q2_num_traits == 2).sum()
            props_q2[i, 2] += (q2_num_traits == 3).sum()
            props_q2[i, 3] += (q2_num_traits > 3).sum()

            specificities[i, 1] += specificity[j, q2_keep].mean() / num_traits
            frequencies[i, 1] += mafs[q2_keep].mean() / num_traits

            q3_keep = (pvals[j] >= b2)*(pvals[j] < b3)
            q3_num_traits = (pvals[:, q3_keep] < sig_p).sum(axis=0)
            this_q3.append((q3_num_traits-1).mean())
            props_q3[i, 0] += (q3_num_traits == 1).sum()
            props_q3[i, 1] += (q3_num_traits == 2).sum()
            props_q3[i, 2] += (q3_num_traits == 3).sum()
            props_q3[i, 3] += (q3_num_traits > 3).sum()

            specificities[i, 2] += specificity[j, q3_keep].mean() / num_traits
            frequencies[i, 2] += mafs[q3_keep].mean() / num_traits

            q4_keep = (pvals[j] >= b3)*(pvals[j] < b4)
            q4_num_traits = (pvals[:, q4_keep] < sig_p).sum(axis=0)
            this_q4.append((q4_num_traits-1).mean())
            props_q4[i, 0] += (q4_num_traits == 1).sum()
            props_q4[i, 1] += (q4_num_traits == 2).sum()
            props_q4[i, 2] += (q4_num_traits == 3).sum()
            props_q4[i, 3] += (q4_num_traits > 3).sum()

            specificities[i, 3] += specificity[j, q4_keep].mean() / num_traits
            frequencies[i, 3] += mafs[q4_keep].mean() / num_traits

        q1.append(np.mean(this_q1))
        q2.append(np.mean(this_q2))
        q3.append(np.mean(this_q3))
        q4.append(np.mean(this_q4))

    df = pd.DataFrame()
    df['rep'] = np.arange(len(props_q1))
    df['p_bin_1_num_hits_1_trait'] = props_q1[:, 0]
    df['p_bin_1_num_hits_2_trait'] = props_q1[:, 1]
    df['p_bin_1_num_hits_3_trait'] = props_q1[:, 2]
    df['p_bin_1_num_hits_4plus_trait'] = props_q1[:, 3]
    df['p_bin_2_num_hits_1_trait'] = props_q2[:, 0]
    df['p_bin_2_num_hits_2_trait'] = props_q2[:, 1]
    df['p_bin_2_num_hits_3_trait'] = props_q2[:, 2]
    df['p_bin_2_num_hits_4plus_trait'] = props_q2[:, 3]
    df['p_bin_3_num_hits_1_trait'] = props_q3[:, 0]
    df['p_bin_3_num_hits_2_trait'] = props_q3[:, 1]
    df['p_bin_3_num_hits_3_trait'] = props_q3[:, 2]
    df['p_bin_3_num_hits_4plus_trait'] = props_q3[:, 3]
    df['p_bin_4_num_hits_1_trait'] = props_q4[:, 0]
    df['p_bin_4_num_hits_2_trait'] = props_q4[:, 1]
    df['p_bin_4_num_hits_3_trait'] = props_q4[:, 2]
    df['p_bin_4_num_hits_4plus_trait'] = props_q4[:, 3]
    df['p_bin_1_mean_num_hits'] = np.array(q1)+1
    df['p_bin_2_mean_num_hits'] = np.array(q2)+1
    df['p_bin_3_mean_num_hits'] = np.array(q3)+1
    df['p_bin_4_mean_num_hits'] = np.array(q4)+1

    df['p_bin_1_specificity'] = specificities[:, 0]
    df['p_bin_2_specificity'] = specificities[:, 1]
    df['p_bin_3_specificity'] = specificities[:, 2]
    df['p_bin_4_specificity'] = specificities[:, 3]

    df['p_bin_1_freq'] = frequencies[:, 0]
    df['p_bin_2_freq'] = frequencies[:, 1]
    df['p_bin_3_freq'] = frequencies[:, 2]
    df['p_bin_4_freq'] = frequencies[:, 3]
    return df


N = 10000000
p = 0.5
sig_p = 1e-5
scale_factors = 0.33

print('\t Default params')
fname = 'data/p_val_ranking_sims.csv'
if not os.path.isfile(fname):
    df = run_sims(N, p, sig_p, scale_factors)
    df.to_csv(fname)

skellington = 'data/p_val_ranking_sims_N_{}_p_{}_sig_p_{}_scale_factor_{}.csv'
for N_var in [N/2, N*2, N*10]:
    print('\tN', N_var)
    fname = skellington.format(
        N_var, p, sig_p, scale_factors
    )
    if os.path.isfile(fname):
        continue
    df = run_sims(N_var, p, sig_p, scale_factors)
    df.to_csv(fname, index=False)

for p_var in [0.1, 0.25, 0.75, 1. + 1./20.]:
    print('\tp', p_var)
    fname = skellington.format(
        N, p_var, sig_p, scale_factors
    )
    if os.path.isfile(fname):
        continue
    df = run_sims(N, p_var, sig_p, scale_factors)
    df.to_csv(fname , index=False)

for sig_p_var in [1e-3, 1e-4, 1e-6, 1e-7]:
    print('\tsig_p', sig_p_var)
    fname = skellington.format(
        N, p, sig_p_var, scale_factors
    )
    if os.path.isfile(fname):
        continue
    df = run_sims(N, p, sig_p_var, scale_factors)
    df.to_csv(fname, index=False)

for scale_factors_var in [
    scale_factors/100, scale_factors/10, scale_factors*10, scale_factors*100
]:
    print('\tscale_factors', scale_factors_var)
    fname = skellington.format(
        N, p, sig_p, scale_factors_var
    )
    if os.path.isfile(fname):
        continue
    df = run_sims(N, p, sig_p, scale_factors_var)
    df.to_csv(fname, index=False)
