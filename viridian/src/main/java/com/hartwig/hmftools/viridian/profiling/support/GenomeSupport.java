package com.hartwig.hmftools.viridian.profiling.support;

import com.hartwig.hmftools.viridian.profiling.common.GenomeStats;
import com.hartwig.hmftools.viridian.profiling.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.jetbrains.annotations.Nullable;

// Per-virus-genome support information from aligning reads to every virus genome at once (BWA-MEM -a mode).
// Used to find which viruses are present in reads and to select the representative for a virus group.
public record GenomeSupport(
        GenomeStats stats,
        // Reads with more than one alignment to this virus genome.
        int multiAlignReads,
        // Alignments to this virus genome per read. Theoretically can be > 1. Null only when every alignment was
        // origin clipped.
        @Nullable SummaryStats alignmentsPerRead,
        // Read attribution taking into account alignment edit distance.
        // Each read gives total 1 vote, spread across its aligned genomes by edit distance.
        // A genome with many votes has many aligned reads.
        double readVotes,
        // Reads whose best alignment, among all genomes, is to this genome.
        // A read counts towards multiple genomes if its alignments are tied.
        int explainedReads)
{
    public VirusGenome genome()
    {
        return stats.genome();
    }
}
