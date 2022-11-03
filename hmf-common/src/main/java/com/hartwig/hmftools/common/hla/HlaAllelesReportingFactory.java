package com.hartwig.hmftools.common.hla;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class HlaAllelesReportingFactory {

    private static final Logger LOGGER = LogManager.getLogger(HlaAllelesReportingFactory.class);

    public static final DecimalFormat SINGLE_DIGIT = new DecimalFormat("#.#", DecimalFormatSymbols.getInstance(Locale.ENGLISH));

    @NotNull
    public static Map<String, List<LilacAllele>> generateLilacMap(@NotNull LilacSummaryData lilacSummaryData) {
        Map<String, List<LilacAllele>> mapLilacReportingAlleles = Maps.newHashMap();
        for (LilacAllele lilacAllele : lilacSummaryData.alleles()) {

            if (mapLilacReportingAlleles.containsKey(lilacAllele.allele())) {
                List<LilacAllele> curent = mapLilacReportingAlleles.get(lilacAllele.allele());
                curent.add(lilacAllele);
                mapLilacReportingAlleles.put(lilacAllele.allele(), curent);
            } else {
                mapLilacReportingAlleles.put(lilacAllele.allele(), Lists.newArrayList(lilacAllele));
            }
        }
        return mapLilacReportingAlleles;
    }

    @NotNull
    public static HlaAllelesReportingData convertToReportData(@NotNull LilacSummaryData lilacSummaryData, boolean hasRealiablePurity) {
        List<HlaReporting> lilacReportingList = Lists.newArrayList();

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

            lilacReportingList.add(ImmutableHlaReporting.builder()
                    .hlaAllele(ImmutableHlaAllele.builder()
                            .germlineAllele(keyMap.getValue().get(0).allele())
                            .gene(extractHLAGene(keyMap.getValue().get(0)))
                            .build())
                    .germlineCopies(germlineCopies)
                    .tumorCopies(tumorCopies)
                    .somaticMutations(mutationString)
                    .interpretation(HLApresenceInTumor(keyMap.getValue().get(0), mutationString, hasRealiablePurity))
                    .build());
        }

        Map<String, List<HlaReporting>> hlaAlleleMap = Maps.newHashMap();

        for (HlaReporting hlaReporting : lilacReportingList) {
            if (hlaAlleleMap.containsKey(hlaReporting.hlaAllele().gene())) {
                List<HlaReporting> curent = hlaAlleleMap.get(hlaReporting.hlaAllele().gene());
                curent.add(hlaReporting);
                hlaAlleleMap.put(hlaReporting.hlaAllele().gene(), curent);
            } else {
                hlaAlleleMap.put(hlaReporting.hlaAllele().gene(), Lists.newArrayList(hlaReporting));
            }
        }

        return ImmutableHlaAllelesReportingData.builder()
                .hlaQC(lilacSummaryData.qc())
                .hlaAllelesReporting(hlaAlleleMap)
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

        if (Doubles.positive(allele.somaticInframeIndel())) {
            joiner.add(SINGLE_DIGIT.format(allele.somaticInframeIndel()) + " inframe indel");
        }

        String result = joiner.toString();
        return !result.isEmpty() ? result : "No";
    }

    @NotNull
    public static String HLApresenceInTumor(@NotNull LilacAllele allele, @NotNull String mutationString, boolean hasReliablePurity) {
        double tumorCopies = Double.parseDouble(SINGLE_DIGIT.format(allele.tumorCopyNumber()));
        boolean mutation = mutationString.contains("missense") || mutationString.contains("nonsense or frameshift")
                || mutationString.contains("splice") || mutationString.contains("inframe indel");
        if (hasReliablePurity) {
            if (tumorCopies == 0) {
                if (mutation) {
                    return "Unknown";
                } else {
                    return "No";
                }
            } else {
                if (!mutation) {
                    return "Yes";
                } else {
                    return "Yes, but mutation(s) detected";
                }
            }
        } else {
            return "Unknown";
        }
    }
}