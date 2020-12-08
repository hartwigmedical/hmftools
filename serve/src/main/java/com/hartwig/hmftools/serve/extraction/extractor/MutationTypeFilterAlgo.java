package com.hartwig.hmftools.serve.extraction.extractor;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class MutationTypeFilterAlgo {

    private static final Logger LOGGER = LogManager.getLogger(MutationTypeFilterAlgo.class);

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
            return MutationTypeFilter.MISSENSE_INFRAME_ANY;
        } else if (formattedEvent.contains("deletion")) {
            return MutationTypeFilter.MISSENSE_INFRAME_DELETION;
        } else if (formattedEvent.contains("insertion")) {
            return MutationTypeFilter.MISSENSE_INFRAME_INSERTION;
        } else if (formattedEvent.contains("frameshift")) {
            return MutationTypeFilter.NONSENSE_OR_FRAMESHIFT;
        } else {
            DriverGene driverGene = findByGene(gene);
            if (driverGene != null) {
                if (driverGene.likelihoodType() == DriverCategory.ONCO) {
                    return MutationTypeFilter.MISSENSE_ANY;
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
