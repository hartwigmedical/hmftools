package com.hartwig.hmftools.viridian.detection.support;

import com.hartwig.hmftools.viridian.detection.common.GenomeStats;
import com.hartwig.hmftools.viridian.detection.common.SummaryStats;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

import org.jetbrains.annotations.Nullable;

// Per-contig support information from aligning reads to every virus genome at once (BWA-MEM -a mode).
// Used to detect viruses and select the representative contig for an oncology group.
public record GenomeSupport(
        GenomeStats stats,
        GenomeFilterStatus filterStatus,
        // Reads with more than one alignment to this contig.
        int multiAlignReads,
        // Alignments to this contig per read. Theoretically can be > 1. Null only when every alignment was origin clipped.
        @Nullable SummaryStats alignmentsPerRead,
        // Read attribution taking into account alignment edit distance (divergence).
        double readVotes)
{
    public VirusGenome genome()
    {
        return stats.genome();
    }

    public boolean isCandidate()
    {
        return filterStatus == GenomeFilterStatus.CANDIDATE;
    }
}
