import os
from .load import *
import seaborn as sns
import matplotlib.pyplot as plt
from .timer import Timer


def plot_all(samples, output_prefix='figures/', fuzzy_distance=15):

    # TODO: add config parsing for plot options

    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    mapped_counts = prepare_mapped_counts(samples)

    df_ref = load_reference_dna()
    unique_names = df_ref['name'].drop_duplicates(keep=False)

    source_flag = reference_contains_source()

    if source_flag:
        
        for field, matches in mapped_counts.items():
            plot_crossmapping(field, matches, output_prefix, 
                unique_names)

        
            plot_abundances_by_sample(field, matches, 
                output_prefix, df_ref)


def divide_matches(matches, unique_names):
    # higher is better
    priority = (matches
    .groupby(['sample', 'source'])
    ['read_count'].sum().rename('priority')
    .reset_index()
    )

    # keep matches to the dominant source 
    # (removes matches to duplicates across sources)
    dominant_matches = (matches
    .merge(priority)
    .sort_values('priority', ascending=False)
    .drop_duplicates(['sample', 'name'])
    )

    unique_matches = (matches
     .loc[lambda x: x['name'].isin(unique_names)]
     .assign(read_fraction=lambda x: 
        x['read_fraction'] 
        / x.groupby('sample')['read_fraction'].transform('sum'))
    )

    return dominant_matches, unique_matches


def plot_crossmapping(field, matches, output_prefix, unique_names):
    dominant_matches, unique_matches = divide_matches(
        matches, unique_names)

    msg = f'Plotting sample vs. source for {field}...'
    with Timer(verbose=msg):
        fig, ax, df_plot = heatmap_sample_vs_source(dominant_matches)
        ax.set_title(f'Read fraction for\n{field}')
        f = f'{output_prefix}heatmap_sample_vs_source_{field}'
        fig.savefig(f, bbox_inches='tight')
        plt.close(fig)
        df_plot.to_csv(f + '.csv', index=None)

        fig, ax, df_plot = heatmap_sample_vs_source(unique_matches)
        ax.set_title(f'Unique read fraction for\n{field}')
        f = f'{output_prefix}heatmap_sample_vs_source_{field}_unique'
        fig.savefig(f, bbox_inches='tight')
        plt.close(fig)
        df_plot.to_csv(f + '.csv', index=None)


def plot_abundances_by_sample(field, matches, output_prefix, 
        df_reference):
    msg = f'Plotting abundance by sample for {field}...'
    with Timer(verbose=msg):
        source_counts = (df_reference
            .drop_duplicates(['source', 'name'])
            .groupby('source').size())
        fig = facet_abundances_by_sample(matches, source_counts)
        f = f'{output_prefix}abundance_by_sample_{field}'
        fig.savefig(f, bbox_inches='tight')
        plt.close(fig)
        matches.to_csv(f + '.csv', index=None)


def prepare_mapped_counts(samples, minimum_distance=None):
    df_mapped, fields = load_data_for_plotting(samples, simulate=True)

    if minimum_distance is None:
        minimum_distance = {}

    df_reference = load_reference_dna().drop('reference_dna', axis=1)

    cols = ['sample', 'name']
    if reference_contains_source():
        cols += ['source']
        
    # prepare data
    mapped_counts = {}
    for field in fields:
        cutoff = minimum_distance.get(field, 0)
        keep = df_mapped[f'{field}_distance'] <= cutoff
        mapped_counts[field] = (df_mapped[keep]
        .rename(columns={f'{field}_match': 'name'})
        .groupby(['sample', 'name'])
        .size().rename('read_count').reset_index()
        .assign(read_fraction=lambda x: 
                x['read_count'] / x.groupby('sample')['read_count'].transform('sum'))
        .merge(df_reference)
        )
    
    return mapped_counts


def heatmap_sample_vs_source(matches):
    df_plot = (matches
     .pivot_table(index='sample', columns='source', 
                  values='read_fraction', aggfunc='sum')
     .fillna(0))
    
    padding = np.array([1, 1.5])
    scaling = np.array([0.5, 0.5])
    figsize = padding + scaling * df_plot.shape
    fig, ax = plt.subplots(figsize=figsize)
    (df_plot
     .pipe(sns.heatmap, annot=True, xticklabels=True, yticklabels=True, 
           cbar=False, ax=ax)
    )
    plt.yticks(rotation=0)
    plt.xticks(rotation=30)
    
    return fig, ax, df_plot


def facet_abundances_by_sample(matches, source_counts):

    def plot(data, label, color):
        data = sorted(data['read_count'])[::-1]
        ax = plt.gca()
        ax.plot(np.arange(len(data)) + 1, data, color=color, label=label, marker='.')
        x = source_counts[label]
        ax.plot([x, x], [1, 100], color=color, ls=':', label='# designs')

    hue_order = natsorted(set(matches['source']))
    row_order = natsorted(set(matches['sample']))
    fg = (matches
    .pipe(sns.FacetGrid, 
          row='sample', row_order=row_order, 
          hue='source', hue_order=hue_order,
          height=2, aspect=2, sharex=False)
    .map_dataframe(plot)
    .add_legend()
    )

    fg.axes.flat[-1].set_xlabel('Rank');
    for ax in fg.axes.flat[:]:
        ax.set_ylabel('Read count')
        ax.set_yscale('log')
        
    return fg.fig