package com.hartwig.hmftools.datamodel.finding;

import java.util.List;
import java.util.Set;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class GeneListUtil
{

    static Set<String> genes(@NotNull List<PurpleVariant> reportableVariants,
            @NotNull List<PurpleGainDeletion> gainDeletions,
            @Nullable List<LinxHomozygousDisruption> homozygousDisruptions,
            @NotNull Set<String> genes)
    {
        Set<String> genesDisplay = Sets.newTreeSet();

        genesDisplay.addAll(filteredMapped(reportableVariants,
                variant -> genes.contains(variant.gene()),
                PurpleVariant::gene));

        genesDisplay.addAll(filteredMapped(gainDeletions,
                gainDeletion -> genes.contains(gainDeletion.gene()) && (
                        gainDeletion.interpretation() == CopyNumberInterpretation.PARTIAL_DEL
                                || gainDeletion.interpretation() == CopyNumberInterpretation.FULL_DEL),
                PurpleGainDeletion::gene));

        if(homozygousDisruptions != null)
        {
            genesDisplay.addAll(filteredMapped(homozygousDisruptions,
                    homozygousDisruption -> genes.contains(homozygousDisruption.gene()),
                    LinxHomozygousDisruption::gene));
        }

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
