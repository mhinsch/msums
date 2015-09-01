#include "runanalysis.h"


void analyse(ostream & abcstatfile, int dataset, 
	vector<vector<SingleStats::sample_t> > & sequences,
	const vector<bool> & show_pops, const vector<vector<bool> > & show_pairs,
	const vector<SingleStats*> & single_stats,
	const vector<PairStats*> & pair_stats,
	const vector<GroupStats*> & group_stats,
	const vector<vector<size_t> > & groups)
	{
	const size_t n_loci = sequences[0].size();
	const size_t n_pops = sequences.size();

	for(size_t l=0; l<n_loci; l++)	
		{
		// new locus
		abcstatfile << "\n" << dataset << "\t" << l;

		for (size_t p=0; p<n_pops; p++)
			{
			if (!show_pops[p])
				continue;

			sequences[p][l].prepare_alleles_per_site();

			for (size_t stat=0; stat<single_stats.size(); stat++)
				abcstatfile << '\t' <<
					single_stats[stat]->analyse(sequences[p][l]);
			}
		
		for (size_t pop1=0; pop1<n_pops-1; pop1++)
			for (size_t pop2=pop1+1; pop2<n_pops; pop2++)
				{
				if (!show_pairs[pop1][pop2])
					continue;

				sequences[pop1][l].prepare_alleles_per_site();
				sequences[pop2][l].prepare_alleles_per_site();

				const SingleStats::sample_t p1 = sequences[pop1][l];
				const SingleStats::sample_t p2 = sequences[pop2][l];

				const PairStats::sample_t ps(p1, p2);
				for (size_t stat=0; stat<pair_stats.size(); stat++)
					abcstatfile << '\t' << pair_stats[stat]->analyse(ps);
				}

		GroupStats::sample_t group;
		for (size_t stat=0; stat<group_stats.size(); stat++)
			{
			group.clear();
			for (size_t i=0; i<groups[stat].size(); i++)
				group.push_back(&(sequences[groups[stat][i]][l]));

			abcstatfile << '\t' << group_stats[stat]->analyse(group);
			}
		}
	}


void analyse_aggr(ostream & abcstatfile, int dataset,
	vector<vector<SingleStats::sample_t> > & sequences,
	const vector<bool> & show_pops, const vector<vector<bool> > & show_pairs,
	const vector<SingleStats*> & single_stats,
	const vector<PairStats*> & pair_stats,
	const vector<GroupStats*> & group_stats,
	const vector<vector<size_t> > & groups)
	{
	const size_t n_pops = sequences.size();

	abcstatfile << dataset;

	for (size_t p=0; p<n_pops; p++)
		{
		if (!show_pops[p]) continue;

		for (size_t stat=0; stat<single_stats.size(); stat++)
			{
			Aggregate aggregate;
			for (size_t l=0; l<sequences[p].size(); l++)
				{
				sequences[p][l].prepare_alleles_per_site();
				aggregate(single_stats[stat]->analyse(sequences[p][l]));
				}

			aggregate.analyse();
			abcstatfile << '\t' << aggregate.mean() << '\t'<< aggregate.std();
			}
		}

	vector<Aggregate> aggregates;

	for (size_t pop1=0; pop1<n_pops-1; pop1++)
		for (size_t pop2=pop1+1; pop2<n_pops; pop2++)
			{
			if (!show_pairs[pop1][pop2])
				continue;

			vector<SingleStats::sample_t> & p1 = sequences[pop1];
			vector<SingleStats::sample_t> & p2 = sequences[pop2];

			aggregates.clear();
			aggregates.resize(pair_stats.size());

			for (size_t l=0; l<p1.size(); l++)
				{
				p1[l].prepare_alleles_per_site();
				p2[l].prepare_alleles_per_site();
				// PairSample caches analysis results
				const PairStats::sample_t ps(p1[l], p2[l]);
				for (size_t stat=0; stat<pair_stats.size(); stat++)
					aggregates[stat](pair_stats[stat]->analyse(ps));
				}

			for (size_t stat=0; stat<pair_stats.size(); stat++)
				{
				aggregates[stat].analyse();
				abcstatfile << '\t' << aggregates[stat].mean() <<
				   	'\t' << aggregates[stat].std();
				}
			}

	GroupStats::sample_t group;
	for (size_t stat=0; stat<group_stats.size(); stat++)
		{
		Aggregate aggregate;
		const vector<size_t> cur_group = groups[stat];
		const GroupStats & ana = *group_stats[stat];

		for (size_t l=0; l<sequences[0].size(); l++)
			{
			group.clear();
			for (size_t i=0; i<groups[stat].size(); i++)
				group.push_back(&(sequences[cur_group[i]][l]));

			aggregate(ana.analyse(group));
			}

		aggregate.analyse();
		abcstatfile << '\t' << aggregate.mean() << '\t'<< aggregate.std();
		}
	}
