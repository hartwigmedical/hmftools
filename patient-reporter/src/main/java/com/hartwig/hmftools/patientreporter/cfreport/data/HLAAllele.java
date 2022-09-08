package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.text.DecimalFormat;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.hla.LilacAllele;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.patientreporter.cfreport.ReportResources;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class HLAAllele {

    private HLAAllele() {
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

    @NotNull
    public static String HLApresenceInTumor(@NotNull LilacAllele allele, @NotNull String mutationString, boolean hasReliablePurity) {
        double tumorCopies = Double.parseDouble(HLAAllele.SINGLE_DIGIT.format(allele.tumorCopyNumber()));
        String presenceInTumor = Strings.EMPTY;
        boolean mutation = mutationString.contains("missense") || mutationString.contains("nonsense or frameshift")
                || mutationString.contains("splice") || mutationString.contains("synonymous") || mutationString.contains("inframe indel");
        if (!hasReliablePurity || (tumorCopies == 0 && !mutationString.equals("None"))) {
            presenceInTumor = "Unknown";
        } else if (tumorCopies >= 1 && mutation) {
            presenceInTumor = "Yes, but mutation(s) detected";
        } else if (tumorCopies == 0 && mutation) {
            presenceInTumor = "Yes";
        } else if (mutationString.equals("None")) {
            presenceInTumor = "No";
        }
        return presenceInTumor;
    }
}