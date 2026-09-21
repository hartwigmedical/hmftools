package com.hartwig.hmftools.viridian.detection.contig_support;

import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_CONTIG_COVERAGE_MIN;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_CONTIG_COVERAGE_MIN_LOWER;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_CONTIG_VOTES_PER_BASE_MIN;

import java.util.Map;

import com.hartwig.hmftools.viridian.reference.ViralContig;

// Determines whether a contig carries enough evidence to be considered present in the sample.
public class ContigSupportFilter
{
    public static Map<ViralContig, ContigFilterStatus> statuses(
            Map<ViralContig, Integer> coveredBases, Map<ViralContig, Double> readVotes, double meanReadLength)
    {
        boolean groupPresent = coveredBases.entrySet().stream()
                .anyMatch(entry -> coverageFraction(entry.getKey(), entry.getValue()) >= VIRAL_CONTIG_COVERAGE_MIN);

        return coveredBases.entrySet().stream().collect(toMap(
                Map.Entry::getKey,
                entry -> status(entry.getKey(), entry.getValue(), readVotes.get(entry.getKey()), groupPresent, meanReadLength)));
    }

    private static ContigFilterStatus status(
            ViralContig contig, int coveredBases, double readVotes, boolean groupPresent, double meanReadLength)
    {
        // Once some contig has established the group, its siblings are kept down to a slightly lower coverage, so a
        // near-identical sibling straddling the cutoff is not harshly lost.
        if(!groupPresent || coverageFraction(contig, coveredBases) < VIRAL_CONTIG_COVERAGE_MIN_LOWER)
        {
            return ContigFilterStatus.LOW_COVERAGE;
        }
        else if(readVotes < voteFloor(contig, meanReadLength))
        {
            return ContigFilterStatus.LOW_VOTE_DENSITY;
        }
        else
        {
            return ContigFilterStatus.CANDIDATE;
        }
    }

    private static double coverageFraction(ViralContig contig, int coveredBases)
    {
        return (double) coveredBases / contig.length();
    }

    // Enough votes to have covered the minimum coverage fraction of the contig at the required depth.
    // Drops contigs which are similar to other contigs and so get coverage, but the reads decisively prefer the other
    // contig.
    private static double voteFloor(ViralContig contig, double meanReadLength)
    {
        return VIRAL_CONTIG_VOTES_PER_BASE_MIN * VIRAL_CONTIG_COVERAGE_MIN * contig.length() / meanReadLength;
    }
}
