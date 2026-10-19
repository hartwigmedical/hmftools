package com.hartwig.hmftools.viridian.detection.support;

import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRUS_GENOME_COVERAGE_MIN;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRUS_GENOME_COVERAGE_MIN_LOWER;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRUS_GENOME_VOTES_PER_BASE_MIN;

import java.util.Map;

import com.hartwig.hmftools.viridian.detection.common.GenomeStats;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Determines whether a virus genome has enough support to be a candidate for detection.
// These filters provide a floor on virus detection sensitivity. Genomes filtered here are obviously marginal.
public class GenomeSupportFilter
{
    public static Map<VirusGenome, GenomeFilterStatus> statuses(
            Map<VirusGenome, GenomeStats> genomeStats, Map<VirusGenome, Double> readVotes, double meanReadLength)
    {
        boolean groupReachesCoverageMin = genomeStats.values().stream()
                .anyMatch(stats -> stats.coverageFraction() >= VIRUS_GENOME_COVERAGE_MIN);

        return genomeStats.entrySet().stream().collect(toMap(
                Map.Entry::getKey,
                entry -> status(
                        entry.getValue(), readVotes.getOrDefault(entry.getKey(), 0.0), groupReachesCoverageMin,
                        meanReadLength)));
    }

    private static GenomeFilterStatus status(
            GenomeStats stats, double readVotes, boolean groupReachesCoverageMin, double meanReadLength)
    {
        // Once some virus genome reaches the coverage minimum, its siblings are kept down to a slightly lower
        // coverage, so a near-identical sibling straddling the cutoff is not harshly lost.
        if(!groupReachesCoverageMin || stats.coverageFraction() < VIRUS_GENOME_COVERAGE_MIN_LOWER)
        {
            return GenomeFilterStatus.LOW_COVERAGE;
        }
        else if(readVotes < voteFloor(stats.genome(), meanReadLength))
        {
            return GenomeFilterStatus.LOW_VOTE_DENSITY;
        }
        else
        {
            return GenomeFilterStatus.PASS;
        }
    }

    // Enough votes to have covered the minimum coverage fraction of the virus genome at the required depth.
    // Drops virus genomes which are similar to others and so get coverage, but the reads decisively prefer the
    // other genome.
    private static double voteFloor(VirusGenome genome, double meanReadLength)
    {
        return VIRUS_GENOME_VOTES_PER_BASE_MIN * VIRUS_GENOME_COVERAGE_MIN * genome.length() / meanReadLength;
    }
}
