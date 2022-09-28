package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.hla.LilacReporting;
import com.hartwig.hmftools.common.utils.DataUtil;

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

    @NotNull
    public static String copyNumberStringTumor(Double copyNumber, boolean hasReliablePurity) {
        return hasReliablePurity && !copyNumber.isNaN() ? String.valueOf(Math.round(copyNumber)) : DataUtil.NA_STRING;
    }

    @NotNull
    public static String copyNumberStringGermline(Double copyNumber, boolean hasReliablePurity) {
        return hasReliablePurity && !copyNumber.isNaN() ? String.valueOf(Math.round(copyNumber)) : DataUtil.NA_STRING;
    }
}