package com.hartwig.hmftools.viridian.detection.support;

import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_CONTIG_COVERAGE_MIN;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_CONTIG_COVERAGE_MIN_LOWER;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_CONTIG_VOTES_PER_BASE_MIN;

import java.util.Map;

import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.reference.ViralContig;

// Determines whether a contig carries enough evidence to be a candidate virus genome.
public class ContigSupportFilter
{
    public static Map<ViralContig, ContigFilterStatus> statuses(
            Map<ViralContig, ContigStats> contigStats, Map<ViralContig, Double> readVotes, double meanReadLength)
    {
        boolean groupReachesCoverageMin = contigStats.values().stream()
                .anyMatch(stats -> stats.coverageFraction() >= VIRAL_CONTIG_COVERAGE_MIN);

        return contigStats.entrySet().stream().collect(toMap(
                Map.Entry::getKey,
                entry -> status(
                        entry.getValue(), readVotes.getOrDefault(entry.getKey(), 0.0), groupReachesCoverageMin,
                        meanReadLength)));
    }

    private static ContigFilterStatus status(
            ContigStats stats, double readVotes, boolean groupReachesCoverageMin, double meanReadLength)
    {
        // Once some contig reaches the coverage minimum, its siblings are kept down to a slightly lower coverage, so a
        // near-identical sibling straddling the cutoff is not harshly lost.
        if(!groupReachesCoverageMin || stats.coverageFraction() < VIRAL_CONTIG_COVERAGE_MIN_LOWER)
        {
            return ContigFilterStatus.LOW_COVERAGE;
        }
        else if(readVotes < voteFloor(stats.contig(), meanReadLength))
        {
            return ContigFilterStatus.LOW_VOTE_DENSITY;
        }
        else
        {
            return ContigFilterStatus.CANDIDATE;
        }
    }

    // Enough votes to have covered the minimum coverage fraction of the contig at the required depth.
    // Drops contigs which are similar to other contigs and so get coverage, but the reads decisively prefer the other
    // contig.
    private static double voteFloor(ViralContig contig, double meanReadLength)
    {
        return VIRAL_CONTIG_VOTES_PER_BASE_MIN * VIRAL_CONTIG_COVERAGE_MIN * contig.length() / meanReadLength;
    }
}
