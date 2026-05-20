package com.hartwig.hmftools.finding;

import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.finding.datamodel.Disruption;
import com.hartwig.hmftools.finding.datamodel.GainDeletion;
import com.hartwig.hmftools.finding.datamodel.SmallVariant;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;

class GeneListUtil
{
    static SortedSet<String> genes(DriverFindingList<SmallVariant> smallVariants,
            DriverFindingList<GainDeletion> gainDeletions,
            List<Disruption> germlineHomozygousDisruptions,
            Set<String> genes)
    {
        SortedSet<String> genesDisplay = new TreeSet<>();

        genesDisplay.addAll(filteredMapped(smallVariants.findings(),
                variant -> genes.contains(variant.gene()),
                SmallVariant::gene));

        genesDisplay.addAll(filteredMapped(gainDeletions.findings(),
                gainDeletion -> genes.contains(gainDeletion.gene()) && gainDeletion.isDeletion() && !gainDeletion.isLossOfHeterozygosity(),
                GainDeletion::gene));

        genesDisplay.addAll(filteredMapped(germlineHomozygousDisruptions,
                homozygousDisruption -> genes.contains(homozygousDisruption.gene()),
                Disruption::gene));

        return genesDisplay;
    }

    private static <T> Set<String> filteredMapped(List<T> items, Predicate<T> filter, Function<T, String> mapper)
    {
        return items.stream()
                .filter(filter)
                .map(mapper)
                .collect(Collectors.toSet());
    }
}
