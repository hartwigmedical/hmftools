package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.hla.LilacReporting;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class HLAAllele {

    private HLAAllele() {
    }


//    @NotNull
//    public static List<LilacAllele> sort(@NotNull List<LilacReporting> alleles) {
//        return alleles.stream()
//                .sorted(Comparator.comparing(LilacAllele::allele).thenComparingInt(LilacAllele::refFragments))
//                .collect(Collectors.toList());
//    }


}