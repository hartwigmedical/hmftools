package com.hartwig.hmftools.common.hla;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class LilacReportingFactory {

    private static final Logger LOGGER = LogManager.getLogger(LilacReportingFactory.class);

    public static final DecimalFormat SINGLE_DIGIT = new DecimalFormat("#.#", DecimalFormatSymbols.getInstance(Locale.ENGLISH));

    @NotNull
    public static Map<String, List<LilacAllele>> generateLilacMap(@NotNull LilacSummaryData lilacSummaryData) {
        Map<String, List<LilacAllele>> mapLilacReportingAlleles = Maps.newHashMap();
        for (LilacAllele lilacAllele : lilacSummaryData.alleles()) {
            List<LilacAllele> variantKeys = Lists.newArrayList();

            if (mapLilacReportingAlleles.containsKey(lilacAllele.allele())) {
                variantKeys.addAll(mapLilacReportingAlleles.get(lilacAllele.allele()));
                variantKeys.add(lilacAllele);
                mapLilacReportingAlleles.put(lilacAllele.allele(), variantKeys);
            } else {
                variantKeys.add(lilacAllele);
                mapLilacReportingAlleles.put(lilacAllele.allele(), variantKeys);
            }
        }
        return mapLilacReportingAlleles;
    }

    @NotNull
    public static LilacReportingData convertToReportData(@NotNull LilacSummaryData lilacSummaryData, boolean hasRealiablePurity) {
        List<LilacReporting> lilacReportingList = Lists.newArrayList();

        Map<String, List<LilacAllele>> mapLilacReportingAlleles = generateLilacMap(lilacSummaryData);

        for (Map.Entry<String, List<LilacAllele>> keyMap : mapLilacReportingAlleles.entrySet()) {
            double germlineCopies = 0;
            double tumorCopies = 0;
            String mutationString = Strings.EMPTY;
            if (keyMap.getValue().size() == 1) {
                germlineCopies = 1;
                tumorCopies = keyMap.getValue().get(0).tumorCopyNumber();
                mutationString = mutationString(keyMap.getValue().get(0));
            } else if (keyMap.getValue().size() == 2) {
                germlineCopies = 2;
                tumorCopies = keyMap.getValue().get(0).tumorCopyNumber() + keyMap.getValue().get(1).tumorCopyNumber();
                //Assume only somatic count is added to one of the 2 alleles
                mutationString = mutationString(keyMap.getValue().get(0));
            } else {
                LOGGER.warn("To many hla alleles of allele '{}'", keyMap.getKey());
            }

            lilacReportingList.add(ImmutableLilacReporting.builder()
                    .lilacGermlineAllele(ImmutableLilacGermlineAllele.builder()
                            .germlineAllele(keyMap.getValue().get(0).allele())
                            .gene(extractHLAGene(keyMap.getValue().get(0)))
                            .build())
                    .germlineCopies(germlineCopies)
                    .tumorCopies(tumorCopies)
                    .somaticMutations(mutationString)
                    .interpretation(HLApresenceInTumor(keyMap.getValue().get(0), mutationString, hasRealiablePurity))
                    .build());
        }

        return ImmutableLilacReportingData.builder()
                .lilacQc(lilacSummaryData.qc())
                .lilacReporting(lilacReportingList)
                .build();
    }

    @NotNull
    public static String extractHLAGene(@NotNull LilacAllele lilacAllele) {
        if (lilacAllele.allele().startsWith("A*")) {
            return "HLA-A";
        } else if (lilacAllele.allele().startsWith("B*")) {
            return "HLA-B";
        } else if (lilacAllele.allele().startsWith("C*")) {
            return "HLA-C";
        } else {
            LOGGER.warn("Unknown HLA gene name '{}' present! ", lilacAllele.allele());
            return Strings.EMPTY;
        }
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
        double tumorCopies = Double.parseDouble(SINGLE_DIGIT.format(allele.tumorCopyNumber()));
        String presenceInTumor = Strings.EMPTY;
        boolean mutation = mutationString.contains("missense") || mutationString.contains("nonsense or frameshift")
                || mutationString.contains("splice") || mutationString.contains("synonymous") || mutationString.contains("inframe indel");
        if (hasReliablePurity) {
            if (tumorCopies == 0) {
                if (mutation) {
                    presenceInTumor = "Unknown";
                } else {
                    presenceInTumor = "No";
                }
            } else {
                if (!mutation) {
                    presenceInTumor = "Yes";
                } else {
                    presenceInTumor = "Yes, but mutation(s) detected";
                }
            }
        } else {
            presenceInTumor = "Unknown";
        }
        return presenceInTumor;
    }
}