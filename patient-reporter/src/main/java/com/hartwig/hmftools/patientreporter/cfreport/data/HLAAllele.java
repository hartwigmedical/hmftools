package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;

import org.jetbrains.annotations.NotNull;

public final class HLAAllele {

    private HLAAllele(){
    }

    public static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#.#");

    @NotNull
    public static List<LilacAllele> sort(@NotNull List<LilacAllele> alleles) {
        return alleles.stream()
                .sorted(Comparator.comparing(LilacAllele::allele).thenComparingInt(LilacAllele::refFragments))
                .collect(Collectors.toList());
    }

    @NotNull
    public static String mutationString(@NotNull LilacAllele allele) {
        StringJoiner joiner = new StringJoiner(", ");
        if (Doubles.positive(allele.somaticMissense())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticMissense()) + " missense");
        }

        if (Doubles.positive(allele.somaticNonsenseOrFrameshift())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticNonsenseOrFrameshift()) + " nonsense or frameshift");
        }

        if (Doubles.positive(allele.somaticSplice())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticSplice()) + " splice");
        }

        if (Doubles.positive(allele.somaticSynonymous())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticSynonymous()) + " synonymous");
        }

        if (Doubles.positive(allele.somaticInframeIndel())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticInframeIndel()) + " inframe indel");
        }

        String result = joiner.toString();
        return !result.isEmpty() ? result : "None";
    }
}