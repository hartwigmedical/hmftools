package com.hartwig.hmftools.viridian.profiling.select;

import static java.util.Comparator.comparing;
import static java.util.Comparator.comparingDouble;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.minBy;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.viridian.profiling.support.GenomeSupport;
import com.hartwig.hmftools.viridian.reference.VirusGenome;

// Picks one representative virus genome from each oncology group for downstream calculations.
// Based on observed data, significant novel virus strain and multiple virus occurrence is rare/nonexistent, so this
// selection step is a safe simplification to allow computing coverage/depth without de novo genome assembly.
public class RepresentativeGenomeSelector
{
    public static List<VirusGenome> select(Collection<GenomeSupport> genomeSupports)
    {
        // In observed cases, highest read votes first is an acceptable selector.
        return genomeSupports.stream()
                .filter(support -> support.stats().readCount() > 0)
                .collect(groupingBy(support -> support.genome().oncologyGroup(), minBy(MOST_VOTES_FIRST)))
                .values().stream()
                .map(support -> support.orElseThrow().genome())
                .sorted(comparing(VirusGenome::oncologyGroup))
                .toList();
    }

    private static final Comparator<GenomeSupport> MOST_VOTES_FIRST =
            comparingDouble(GenomeSupport::readVotes).reversed().thenComparing(GenomeSupport::genome);
}
