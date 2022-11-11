package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.hla.HlaReporting;
import com.hartwig.hmftools.common.utils.DataUtil;

import org.jetbrains.annotations.NotNull;

public final class HLAAllele {

    private HLAAllele() {
    }

    @NotNull
    public static List<HlaReporting> sort(@NotNull List<HlaReporting> alleles) {
        return alleles.stream()
                .sorted(Comparator.comparing((HlaReporting lilacReporting) -> lilacReporting.hlaAllele().gene())
                        .thenComparing(germlineAllele -> germlineAllele.hlaAllele().germlineAllele()))
                .collect(Collectors.toList());
    }

    @NotNull
    public static String copyNumberStringTumor(Double copyNumber, boolean hasReliablePurity) {
        return hasReliablePurity && !copyNumber.isNaN() ? String.valueOf(Math.round(copyNumber)) : DataUtil.NA_STRING;
    }

    @NotNull
    public static String copyNumberStringGermline(Double copyNumber, boolean hasReliablePurity) {
        return hasReliablePurity && !copyNumber.isNaN() ? String.valueOf(Math.round(copyNumber)) : DataUtil.NA_STRING;
    }
}