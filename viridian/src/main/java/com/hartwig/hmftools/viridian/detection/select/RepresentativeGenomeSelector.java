package com.hartwig.hmftools.viridian.detection.select;

import static java.util.Comparator.comparing;
import static java.util.Comparator.comparingDouble;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.minBy;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.viridian.detection.support.GenomeSupport;

// Picks one representative virus genome from each oncology group for downstream calculations.
public class RepresentativeGenomeSelector
{
    public static List<RepresentativeGenome> select(Collection<GenomeSupport> genomeSupports)
    {
        return genomeSupports.stream()
                .filter(support -> !support.isFiltered())
                .collect(groupingBy(support -> support.genome().oncologyGroup(), minBy(MOST_VOTES_FIRST)))
                .values().stream()
                .map(support -> new RepresentativeGenome(support.orElseThrow().genome(), RepresentativeStatus.PASS))
                .sorted(comparing(RepresentativeGenome::oncologyGroup))
                .toList();
    }

    private static final Comparator<GenomeSupport> MOST_VOTES_FIRST =
            comparingDouble(GenomeSupport::readVotes).reversed().thenComparing(GenomeSupport::genome);
}
