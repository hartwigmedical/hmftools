package com.hartwig.hmftools.finding;

import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

import org.jetbrains.annotations.Nullable;

class GeneListUtil
{
    static List<String> genes(List<PurpleVariant> reportableVariants,
            List<PurpleGainDeletion> gainDeletions,
            @Nullable List<LinxBreakend> homozygousDisruptions,
            Set<String> genes)
    {
        Set<String> genesDisplay = Sets.newTreeSet();

        genesDisplay.addAll(filteredMapped(reportableVariants,
                variant -> genes.contains(variant.gene()),
                PurpleVariant::gene));

        genesDisplay.addAll(filteredMapped(gainDeletions,
                gainDeletion -> genes.contains(gainDeletion.gene()) && gainDeletion.driver().type() == PurpleDriverType.DEL,
                PurpleGainDeletion::gene));

        if(homozygousDisruptions != null)
        {
            genesDisplay.addAll(filteredMapped(homozygousDisruptions,
                    homozygousDisruption -> genes.contains(homozygousDisruption.gene()),
                    LinxBreakend::gene));
        }

        return genesDisplay.stream().sorted().toList();
    }

    private static <T> Set<String> filteredMapped(List<T> items, Predicate<T> filter, Function<T, String> mapper)
    {
        return items.stream()
                .filter(filter)
                .map(mapper)
                .collect(Collectors.toSet());
    }
}
