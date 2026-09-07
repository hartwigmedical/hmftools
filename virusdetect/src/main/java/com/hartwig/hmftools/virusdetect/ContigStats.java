package com.hartwig.hmftools.virusdetect;

import java.util.Optional;

// Per-contig support over a set of aligned reads, plus how it fares against the sibling contigs of the same virus.
public record ContigStats(
        String contig,
        int contigLength,
        int readCount,              // reads with any alignment to this contig
        int coveredBases,           // contig positions with at least one aligned base
        // Depth spans the whole contig, so uncovered positions count as depth 0.
        SummaryStats depth,
        double meanAlignerScore,    // mean BWA alignment score across those reads
        // Strain rivalry against the near-identical sibling contigs of the same virus:
        double readVotes,           // reads softly attributed to this strain, split across contigs by divergence
        int readsBestInRivals,      // reads that strictly beat every rival contig (ties credit no contig)
        // Divergence lead over the runner-up on strictly-won reads; empty when there were none.
        Optional<SummaryStats> margins)
{
    public double coverageFraction()
    {
        return contigLength == 0 ? 0 : (double) coveredBases / contigLength;
    }
}
