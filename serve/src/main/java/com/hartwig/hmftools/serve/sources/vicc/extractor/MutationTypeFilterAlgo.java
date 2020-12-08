package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.serve.actionability.range.MutationTypeFilter;

import org.jetbrains.annotations.NotNull;

public class MutationTypeFilterAlgo {

    @NotNull
    private final List<DriverGene> driverGenes;

    public MutationTypeFilterAlgo(@NotNull final List<DriverGene> driverGenes) {
        this.driverGenes = driverGenes;
    }

    @NotNull
    public MutationTypeFilter determine(@NotNull String gene, @NotNull String event) {
        String lowerCaseEvent = event.toLowerCase();
        String extractSpecificInfoOfEvent = lowerCaseEvent.substring(lowerCaseEvent.lastIndexOf(" ") + 1);
        MutationTypeFilter filter;
        if (lowerCaseEvent.contains("skipping mutation") || lowerCaseEvent.contains("splice site insertion")) {
            filter = MutationTypeFilter.SPLICE;
        } else if (extractSpecificInfoOfEvent.equals("deletions") || extractSpecificInfoOfEvent.equals("deletion") || lowerCaseEvent.contains(
                "partial deletion of exons")) {
            filter = MutationTypeFilter.MISSENSE_INFRAME_DELETION;
        } else if (extractSpecificInfoOfEvent.equals("insertions") || extractSpecificInfoOfEvent.equals("insertion")) {
            filter = MutationTypeFilter.MISSENSE_INFRAME_INSERTION;
        } else if (extractSpecificInfoOfEvent.equals("deletion/insertion") || extractSpecificInfoOfEvent.equals("insertions/deletions")) {
            filter = MutationTypeFilter.MISSENSE_INFRAME_ANY;
        } else if (extractSpecificInfoOfEvent.equals("frameshift")) {
            filter = MutationTypeFilter.NONSENSE_OR_FRAMESHIFT;
        } else {
            filter = MutationTypeFilter.UNKNOWN;
        }

        for (DriverGene driverGene : driverGenes) {
            if (driverGene.gene().equals(gene)) {
                if (driverGene.likelihoodType() == DriverCategory.ONCO) {
                    if (filter == MutationTypeFilter.UNKNOWN) {
                        return MutationTypeFilter.MISSENSE_ANY;
                    } else {
                        return filter;
                    }
                } else if (driverGene.likelihoodType() == DriverCategory.TSG) {
                    if (filter == MutationTypeFilter.UNKNOWN) {
                        return MutationTypeFilter.ANY;

                    } else {
                        return filter;
                    }
                }
            }
        }

        return filter;
    }
}
