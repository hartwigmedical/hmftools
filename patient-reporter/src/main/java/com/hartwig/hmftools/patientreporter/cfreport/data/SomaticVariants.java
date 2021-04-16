package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.patientreporter.germline.GermlineCondition;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingEntry;
import com.hartwig.hmftools.patientreporter.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.purple.DriverInterpretation;
import com.hartwig.hmftools.protect.purple.ReportableVariant;
import com.hartwig.hmftools.protect.purple.ReportableVariantSource;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SomaticVariants {

    private SomaticVariants() {
    }

    @NotNull
    public static List<ReportableVariant> sort(@NotNull List<ReportableVariant> variants) {
        return variants.stream().sorted((variant1, variant2) -> {
            if (Math.abs(variant1.driverLikelihood() - variant2.driverLikelihood()) > 0.001) {
                return (variant1.driverLikelihood() - variant2.driverLikelihood()) < 0 ? 1 : -1;
            } else {
                if (variant1.gene().equals(variant2.gene())) {
                    // sort on codon position if gene is the same
                    return extractCodonField(variant1.canonicalHgvsCodingImpact()) - extractCodonField(variant2.canonicalHgvsCodingImpact())
                            < 0 ? -1 : 1;
                } else {
                    return variant1.gene().compareTo(variant2.gene());
                }
            }
        }).collect(Collectors.toList());
    }

    public static boolean hasNotifiableGermlineVariant(@NotNull List<ReportableVariant> variants,
            @NotNull GermlineReportingModel germlineReportingModel, @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        for (ReportableVariant variant : variants) {
            if (notifyAboutVariant(variant, germlineReportingModel, germlineReportingLevel)) {
                return true;
            }
        }

        return false;
    }

    @NotNull
    public static String geneDisplayString(@NotNull ReportableVariant variant, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        if (notifyAboutVariant(variant, germlineReportingModel, germlineReportingLevel)) {
            return variant.gene() + " #";
        } else {
            return variant.gene();
        }
    }

    private static boolean notifyAboutVariant(@NotNull ReportableVariant variant, @NotNull GermlineReportingModel germlineReportingModel,
            @NotNull LimsGermlineReportingLevel germlineReportingLevel) {
        boolean includeVariant = false;
        if (variant.source() == ReportableVariantSource.GERMLINE) {
            GermlineReportingEntry reportingEntry = germlineReportingModel.entryForGene(variant.gene());
            if (reportingEntry != null) {
                if (reportingEntry.notifyClinicalGeneticist() == GermlineCondition.ONLY_GERMLINE_HOM) {
                    String conditionFilter = reportingEntry.conditionFilter();
                    if (conditionFilter != null) {
                        includeVariant = variant.genotypeStatus().simplifiedDisplay().equals(conditionFilter);
                    }

                } else if (reportingEntry.notifyClinicalGeneticist() == GermlineCondition.ONLY_SPECIFIC_VARIANT) {
                    String conditionFilter = reportingEntry.conditionFilter();
                    if (conditionFilter != null) {
                        includeVariant = variant.canonicalHgvsProteinImpact().equals(conditionFilter);
                    }
                } else if (reportingEntry.notifyClinicalGeneticist() == GermlineCondition.ALWAYS) {
                    includeVariant = true;
                }
            }
        }

        return includeVariant && germlineReportingModel.notifyAboutGene(variant.gene(), germlineReportingLevel);
    }

    @VisibleForTesting
    static int extractCodonField(@NotNull String hgvsCoding) {
        StringBuilder codonAppender = new StringBuilder();
        boolean noDigitFound = true;
        // hgvsCoding starts with "c.", we need to skip that...
        int index = 2;
        while (noDigitFound && index < hgvsCoding.length()) {
            if ((Character.toString(hgvsCoding.charAt(index)).equals("-") && index == 2) || Character.isDigit(hgvsCoding.charAt(index))) {
                codonAppender.append(hgvsCoding.charAt(index));
            } else {
                noDigitFound = false;
            }
            index++;
        }
        return Integer.parseInt(codonAppender.toString());
    }

    @NotNull
    public static String tVAFString(@NotNull String tVAF, boolean hasReliablePurity) {
        return hasReliablePurity ? tVAF : DataUtil.NA_STRING;
    }

    @NotNull
    public static String copyNumberString(double copyNumber, boolean hasReliablePurity) {
        return hasReliablePurity ? String.valueOf(Math.round(copyNumber)) : DataUtil.NA_STRING;
    }

    @NotNull
    public static String hotspotString(@NotNull Hotspot hotspot) {
        switch (hotspot) {
            case HOTSPOT:
                return "Yes";
            case NEAR_HOTSPOT:
                return "Near";
            default:
                return Strings.EMPTY;
        }
    }

    @NotNull
    public static String biallelicString(@Nullable Boolean biallelic, boolean hasReliablePurity) {
        if (hasReliablePurity) {
            assert biallelic != null;
            return biallelic ? "Yes" : "No";
        } else {
            return DataUtil.NA_STRING;
        }
    }

    @NotNull
    public static String clonalString(double clonalLikelihood) {
        if (clonalLikelihood > 0.95) {
            return ">95%";
        } else if (clonalLikelihood > 0.9) {
            return "90-95%";
        } else if (clonalLikelihood > 0.8) {
            return "80-90%";
        } else if (clonalLikelihood > 0.7) {
            return "70-80%";
        } else if (clonalLikelihood > 0.6) {
            return "60-70%";
        } else if (clonalLikelihood > 0.5) {
            return "50-60%";
        } else if (clonalLikelihood > 0.4) {
            return "40-50%";
        } else if (clonalLikelihood > 0.3) {
            return "30-40%";
        } else if (clonalLikelihood > 0.2) {
            return "20-30%";
        } else if (clonalLikelihood > 0.1) {
            return "10-20%";
        } else if (clonalLikelihood > 0.05) {
            return "5-10%";
        } else {
            return "<5%";
        }
    }

    @NotNull
    public static Set<String> driverGenesWithVariant(@NotNull List<ReportableVariant> variants) {
        Set<String> genes = Sets.newHashSet();
        for (ReportableVariant variant : variants) {
            if (DriverInterpretation.interpret(variant.driverLikelihood()) == DriverInterpretation.HIGH) {
                genes.add(variant.gene());
            }
        }
        return genes;
    }

    public static int countReportableVariants(@NotNull List<ReportableVariant> variants) {
        return variants.size();
    }
}
