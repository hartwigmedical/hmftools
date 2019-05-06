package com.hartwig.hmftools.patientreporter.cfreport.data;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.variant.Clonality;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.patientreporter.cfreport.MathUtil;
import com.hartwig.hmftools.patientreporter.report.data.SomaticVariantDataSource;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import static com.google.common.base.Strings.repeat;

public final class SomaticVariants {

    private SomaticVariants() {
    }

    @NotNull
    public static List<ReportableVariant> sort(@NotNull List<ReportableVariant> variants) {
        return variants.stream().sorted((variant1, variant2) -> {
            Double variant1DriverLikelihood = variant1.driverLikelihood();
            Double variant2DriverLikelihood = variant2.driverLikelihood();

            // Force any variant outside of driver catalog to the bottom of table.
            double driverLikelihood1 = variant1DriverLikelihood != null ? variant1DriverLikelihood : -1;
            double driverLikelihood2 = variant2DriverLikelihood != null ? variant2DriverLikelihood : -1;
            if (Math.abs(driverLikelihood1 - driverLikelihood2) > 0.001) {
                return (driverLikelihood1 - driverLikelihood2) < 0 ? 1 : -1;
            } else {
                if (variant1.gene().equals(variant2.gene())) {
                    // sort on codon position if gene is the same
                    return extractCodonField(variant2.hgvsCodingImpact()).compareTo(extractCodonField(variant1.hgvsCodingImpact()));
                } else {
                    return variant1.gene().compareTo(variant2.gene());
                }
            }
        }).collect(Collectors.toList());
    }

    public static String getGeneDisplayString(@NotNull ReportableVariant variant) {
        return variant.isDrupActionable() ? variant.gene() + " *" : variant.gene();
    }

    @NotNull
    public static String getPloidyVaf(double adjustedCopyNumber, double minorAllelePloidy, double adjustedVAF,
            boolean hasReliablePurityFit) {
        if (!hasReliablePurityFit) {
            return DataUtil.NAString;
        }

        return descriptiveBAF(adjustedCopyNumber, minorAllelePloidy) + " (" + DataUtil.formatPercentage(MathUtil.mapPercentageClamped(
                adjustedVAF,
                0,
                1)) + ")";
    }

    @NotNull
    @VisibleForTesting
    private static String descriptiveBAF(double adjustedCopyNumber, double minorAllelePloidy) {
        int totalAlleleCount = (int) Math.max(0, Math.round(adjustedCopyNumber));
        int minorAlleleCount = (int) Math.max(0, Math.round(minorAllelePloidy));
        int majorAlleleCount = Math.max(0, totalAlleleCount - minorAlleleCount);

        return formatBAFField("A", Math.max(minorAlleleCount, majorAlleleCount)) + formatBAFField("B",
                Math.min(minorAlleleCount, majorAlleleCount));
    }

    @NotNull
    private static String formatBAFField(@NotNull String allele, int count) {
        return count < 10 ? repeat(allele, count) : allele + "[" + count + "x]";
    }

    @NotNull
    private static String extractCodonField(@NotNull String hgvsCoding) {
        StringBuilder stringAppend = new StringBuilder();
        String codonSplit = hgvsCoding.substring(2);
        String codon = "";
        for (int i = 0; i < codonSplit.length(); i++) {
            if (Character.isDigit(codonSplit.charAt(i))) {
                codon = stringAppend.append(codon).append(codonSplit.charAt(i)).toString();
            } else if (!Character.isDigit(codonSplit.charAt(i))) {
                break;
            }
        }
        return codon;
    }

    @NotNull
    public static String getHotspotString(@NotNull Hotspot hotspot) {
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
    public static String getClonalityString(@NotNull Clonality clonality, boolean hasReliablePurityFit) {
        if (!hasReliablePurityFit) {
            return DataUtil.NAString;
        }

        switch (clonality) {
            case CLONAL:
                return "Clonal";
            case SUBCLONAL:
                return "Subclonal";
            case INCONSISTENT:
                return "Inconsistent";
            default:
                return Strings.EMPTY;
        }
    }

    @NotNull
    public static String getBiallelicString(boolean biallelic, DriverCategory driverCategory, boolean hasReliablePurityFit) {
        if (!hasReliablePurityFit) {
            return DataUtil.NAString;
        }

        if (driverCategory != DriverCategory.ONCO) {
            return biallelic ? "Yes" : "No";
        } else {
            return Strings.EMPTY;
        }

    }

    @NotNull
    public static String getDriverString(@Nullable Double driverLikelihood) {
        if (driverLikelihood == null) {
            return Strings.EMPTY;
        }

        if (driverLikelihood > 0.8) {
            return "High";
        } else if (driverLikelihood > 0.2) {
            return "Medium";
        } else {
            return "Low";
        }

    }

    @NotNull
    public static String[] somaticVariantsWithDriver(@NotNull List<ReportableVariant> variants) {
        final List<String> returnVariants = new ArrayList<>();
        for (final ReportableVariant variant : SomaticVariantDataSource.sort(variants)) {
            if (SomaticVariantDataSource.driverField(variant).equals("High")) {
                returnVariants.add(variant.gene());
            }
        }
        return returnVariants.toArray(new String[0]);
    }

    public static int countSomaticVariants(@NotNull List<ReportableVariant> variants) {
        return variants.size();
    }
}
