package com.hartwig.hmftools.viridian.profiling.select;

import static java.util.Comparator.comparing;
import static java.util.Comparator.comparingDouble;
import static java.util.stream.Collectors.groupingBy;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.SECONDARY_GENOME_COUNT;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.viridian.profiling.support.GenomeSupport;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Picks one representative virus genome from each virus group for downstream calculations.
// Based on observed data, significant novel virus strain and multiple virus occurrence is rare/nonexistent, so this
// selection step is a safe simplification to allow computing coverage/depth without de novo genome assembly.
public class RepresentativeGenomeSelector
{
    public static List<RepresentativeSelection> select(Collection<GenomeSupport> genomeSupports)
    {
        return genomeSupports.stream()
                .filter(support -> support.stats().readCount() > 0)
                .collect(groupingBy(support -> support.genome().virusGroup()))
                .values().stream()
                .map(RepresentativeGenomeSelector::selectGroup)
                .sorted(comparing(selection -> selection.representativeGenome().virusGroup()))
                .toList();
    }

    private static RepresentativeSelection selectGroup(List<GenomeSupport> groupSupports)
    {
        // In observed cases, highest read votes first is an acceptable selector.
        List<VirusGenome> ranked = groupSupports.stream()
                .sorted(MOST_VOTES_FIRST)
                .map(GenomeSupport::genome)
                .limit(1 + SECONDARY_GENOME_COUNT)
                .toList();
        return new RepresentativeSelection(ranked.get(0), ranked.subList(1, ranked.size()));
    }

    private static final Comparator<GenomeSupport> MOST_VOTES_FIRST =
            comparingDouble(GenomeSupport::readVotes).reversed().thenComparing(GenomeSupport::genome);
}
