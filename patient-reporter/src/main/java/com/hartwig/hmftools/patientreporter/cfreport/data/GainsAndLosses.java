package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberInterpretation;
import com.hartwig.hmftools.patientreporter.purple.ReportableGainLoss;

import org.jetbrains.annotations.NotNull;

public final class GainsAndLosses {

    private GainsAndLosses() {
    }

    @NotNull
    public static List<ReportableGainLoss> sort(@NotNull List<ReportableGainLoss> reportableGainsAndLosses) {
        return reportableGainsAndLosses.stream().sorted((gainLoss1, gainLoss2) -> {
            String location1 = GeneUtil.zeroPrefixed(gainLoss1.chromosome() + gainLoss1.chromosomeBand());
            String location2 = GeneUtil.zeroPrefixed(gainLoss2.chromosome() + gainLoss2.chromosomeBand());

            if (location1.equals(location2)) {
                return gainLoss1.gene().compareTo(gainLoss2.gene());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static Set<String> amplifiedGenes(@NotNull List<ReportableGainLoss> reportableGainLosses) {
        final Set<String> genes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainLosses) {
            if (gainLoss.interpretation() == CopyNumberInterpretation.GAIN && !gainLoss.gene().isEmpty()) {
                genes.add(gainLoss.gene());
            }
        }
        return genes;
    }

    @NotNull
    public static Set<String> lostGenes(@NotNull List<ReportableGainLoss> reportableGainLosses) {
        final Set<String> genes = Sets.newHashSet();
        for (ReportableGainLoss gainLoss : reportableGainLosses) {
            if (!gainLoss.gene().isEmpty() && (gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS
                    || gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS)) {
                genes.add(gainLoss.gene());
            }
        }
        return genes;
    }
}
