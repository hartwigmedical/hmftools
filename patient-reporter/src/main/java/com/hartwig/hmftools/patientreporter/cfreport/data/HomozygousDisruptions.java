package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;

import org.jetbrains.annotations.NotNull;

public final class HomozygousDisruptions {

    private HomozygousDisruptions() {
    }

    @NotNull
    public static List<HomozygousDisruption> sort(@NotNull List<HomozygousDisruption> homozygousDisruptions) {
        return homozygousDisruptions.stream().sorted((disruption1, disruption2) -> {
            String location1 = GeneUtil.zeroPrefixed(disruption1.chromosome() + disruption1.chromosomeBand());
            String location2 = GeneUtil.zeroPrefixed(disruption2.chromosome() + disruption2.chromosomeBand());

            if (location1.equals(location2)) {
                return disruption1.gene().compareTo(disruption2.gene());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static Set<String> disruptedGenes(@NotNull List<HomozygousDisruption> homozygousDisruptions) {
        Set<String> genes = Sets.newHashSet();
        for (HomozygousDisruption disruption : homozygousDisruptions) {
            genes.add(disruption.gene());
        }
        return genes;
    }
}
