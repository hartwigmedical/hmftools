package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.linx.ReportableGeneDisruption;

import org.jetbrains.annotations.NotNull;

public final class GeneDisruptions {

    private GeneDisruptions() {
    }

    @NotNull
    public static List<ReportableGeneDisruption> sort(@NotNull List<ReportableGeneDisruption> disruptions) {
        return disruptions.stream().sorted((disruption1, disruption2) -> {
            String locationAndGene1 = GeneUtil.zeroPrefixed(disruption1.location()) + disruption1.gene();
            String locationAndGene2 = GeneUtil.zeroPrefixed(disruption2.location()) + disruption2.gene();

            if (locationAndGene1.equals(locationAndGene2)) {
                return disruption1.firstAffectedExon() - disruption2.firstAffectedExon();
            } else {
                return locationAndGene1.compareTo(locationAndGene2);
            }
        }).collect(Collectors.toList());
    }
}
