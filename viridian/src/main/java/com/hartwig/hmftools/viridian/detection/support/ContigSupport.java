package com.hartwig.hmftools.viridian.detection.support;

import com.hartwig.hmftools.viridian.common.SummaryStats;
import com.hartwig.hmftools.viridian.detection.contig_stats.ContigStats;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.jetbrains.annotations.Nullable;

// Per-contig support information from aligning reads to every viral genome at once (BWA-MEM -a mode).
// Used to determine viral presence and select the representative contig for an oncology group.
public record ContigSupport(
        ContigStats stats,
        ContigFilterStatus filterStatus,
        // Reads with more than one alignment to this contig.
        int multiAlignReads,
        // Alignments to this contig per read. Theoretically can be > 1. Null only when every alignment was origin clipped.
        @Nullable SummaryStats alignmentsPerRead,
        // Read attribution taking into account alignment edit distance (divergence).
        double readVotes)
{
    public ViralContig contig()
    {
        return stats.contig();
    }

    public boolean isCandidate()
    {
        return filterStatus == ContigFilterStatus.CANDIDATE;
    }
}
