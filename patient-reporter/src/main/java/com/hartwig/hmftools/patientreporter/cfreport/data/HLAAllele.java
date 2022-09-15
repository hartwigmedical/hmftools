package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.hla.LilacReporting;

import org.jetbrains.annotations.NotNull;

public final class HLAAllele {

    private HLAAllele() {
    }

    @NotNull
    public static List<LilacReporting> sort(@NotNull List<LilacReporting> alleles) {
        return alleles.stream()
                .sorted(Comparator.comparing((LilacReporting germlineAllele) -> germlineAllele.lilacGermlineAllele().gene())
                        .thenComparing(germlineAllele -> germlineAllele.lilacGermlineAllele().germlineAllele()))
                .collect(Collectors.toList());
    }
}