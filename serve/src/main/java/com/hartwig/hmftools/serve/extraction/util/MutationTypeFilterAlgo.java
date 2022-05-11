package com.hartwig.hmftools.serve.extraction.util;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class MutationTypeFilterAlgo {

    @NotNull
    private final List<DriverGene> driverGenes;

    public MutationTypeFilterAlgo(@NotNull final List<DriverGene> driverGenes) {
        this.driverGenes = driverGenes;
    }

    @NotNull
    public MutationTypeFilter determine(@NotNull String gene, @NotNull String event) {
        String formattedEvent = event.toLowerCase();

        if (formattedEvent.contains("skipping") || formattedEvent.contains("splice")) {
            return MutationTypeFilter.SPLICE;
        } else if (formattedEvent.contains("deletion/insertion") || formattedEvent.contains("insertions/deletions")) {
            return MutationTypeFilter.INFRAME;
        } else if (formattedEvent.contains("deletion") || formattedEvent.contains("del")) {
            return MutationTypeFilter.INFRAME_DELETION;
        } else if (formattedEvent.contains("insertion") || formattedEvent.contains("ins")) {
            return MutationTypeFilter.INFRAME_INSERTION;
        } else if (formattedEvent.contains("frameshift")) {
            return MutationTypeFilter.NONSENSE_OR_FRAMESHIFT;
        } else {
            DriverGene driverGene = findByGene(gene);
            if (driverGene != null) {
                if (driverGene.likelihoodType() == DriverCategory.ONCO) {
                    return MutationTypeFilter.MISSENSE;
                } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                    return MutationTypeFilter.ANY;
                }
            }
        }

        return MutationTypeFilter.ANY;
    }

    @Nullable
    private DriverGene findByGene(@NotNull String gene) {
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(gene)) {
                return driverGene;
            }
        }
        return null;
    }
}
