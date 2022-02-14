package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.utils.DataUtil;
import com.hartwig.hmftools.common.variant.DriverInterpretation;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.ReportableVariant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class SomaticVariants {

    private SomaticVariants() {
    }

    private static final String UPSTREAM_GENE_VARIANT = "upstream_gene_variant";

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

    @NotNull
    public static String determineCanonicalImpact(@NotNull String canonicalHgvsCodingImpact, @NotNull String canonicalEffect) {
        String consequence;
        if (canonicalEffect.equals(UPSTREAM_GENE_VARIANT)) {
            consequence = "upstream";
        } else {
            consequence = canonicalHgvsCodingImpact;
        }
        return consequence;
    }

    public static boolean hasNotifiableGermlineVariant(@NotNull Map<ReportableVariant, Boolean> notifyGermlineStatusPerVariant) {
        for (Boolean notify : notifyGermlineStatusPerVariant.values()) {
            if (notify) {
                return true;
            }
        }

        return false;
    }

    @NotNull
    public static String geneDisplayString(@NotNull ReportableVariant variant, boolean notifyGermline) {
        if (notifyGermline) {
            return variant.gene() + " #";
        } else {
            return variant.gene();
        }
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

    @Nullable
    public static String tVAFString(@Nullable String tVAF, boolean hasReliablePurity, Double totalCopyNumber) {
        if (totalCopyNumber == null) {
            return DataUtil.NA_STRING;
        } else {
            double flooredCopyNumber = Math.max(0, totalCopyNumber);
            long roundedCopyNumber = Math.round(flooredCopyNumber);
            return hasReliablePurity && roundedCopyNumber >= 1 ? tVAF : DataUtil.NA_STRING;
        }
    }

    @NotNull
    public static String copyNumberString(Double copyNumber, boolean hasReliablePurity) {
        return hasReliablePurity && !copyNumber.isNaN() ? String.valueOf(Math.round(copyNumber)) : DataUtil.NA_STRING;
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
        if (hasReliablePurity && biallelic != null) {
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

    @NotNull
    public static Set<String> determineMSIgenes(@NotNull List<ReportableVariant> reportableVariants) {
        Set<String> MSI_genes = Sets.newHashSet("MLH1", "MSH2", "MSH6", "PMS2", "EPCAM");
        Set<String> genesDisplay = Sets.newHashSet();

        for (ReportableVariant variant : reportableVariants) {
            if (MSI_genes.contains(variant.gene())) {
                genesDisplay.add(variant.gene());
            }
        }
        return genesDisplay;
    }

    @NotNull
    public static Set<String> determineHRDgenes(@NotNull List<ReportableVariant> reportableVariants) {
        Set<String> HRD_genes = Sets.newHashSet("BRCA1", "BRCA2", "PALB2", "RAD51B", "RAD51D");
        Set<String> genesDisplay = Sets.newHashSet();

        for (ReportableVariant variant : reportableVariants) {
            if (HRD_genes.contains(variant.gene())) {
                genesDisplay.add(variant.gene());
            }
        }
        return genesDisplay;
    }
}
