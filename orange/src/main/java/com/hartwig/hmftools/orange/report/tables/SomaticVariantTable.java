package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.codon.AminoAcidFunctions;
import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class SomaticVariantTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");
    private static final DecimalFormat PERCENTAGE_FORMAT = ReportResources.decimalFormat("#'%'");

    private SomaticVariantTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, @NotNull List<ReportableVariant> driverVariants) {
        if (driverVariants.isEmpty()) {
            return TableUtil.createEmptyTable(title);
        }

        Table table = TableUtil.createReportContentTable(new float[] { 3, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Variant"), TableUtil.createHeaderCell("CN"), TableUtil.createHeaderCell("MACN"),
                        TableUtil.createHeaderCell("VCN"), TableUtil.createHeaderCell("RNA VAF"), TableUtil.createHeaderCell("Biallelic"),
                        TableUtil.createHeaderCell("Hotspot"), TableUtil.createHeaderCell("DL"), TableUtil.createHeaderCell("CL"),
                        TableUtil.createHeaderCell("Phase") });

        for (ReportableVariant variant : sort(driverVariants)) {
            table.addCell(TableUtil.createContentCell(variantField(variant)));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(variant.totalCopyNumber())));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(variant.minorAlleleCopyNumber())));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(variant.alleleCopyNumber())));
            table.addCell(TableUtil.createContentCell("NA"));
            table.addCell(TableUtil.createContentCell(variant.biallelic() ? "Yes" : "No"));
            table.addCell(TableUtil.createContentCell(hotspotField(variant)));
            table.addCell(TableUtil.createContentCell(driverLikelihoodField(variant)));
            table.addCell(TableUtil.createContentCell(PERCENTAGE_FORMAT.format(variant.clonalLikelihood() * 100)));
            table.addCell(TableUtil.createContentCell(
                    variant.localPhaseSet() != null ? String.valueOf(variant.localPhaseSet()) : Strings.EMPTY));
        }

        return TableUtil.createWrappingReportTable(table, title);
    }

    @NotNull
    private static List<ReportableVariant> sort(@NotNull List<ReportableVariant> variants) {
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

    private static int extractCodonField(@NotNull String hgvsCoding) {
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
    private static String variantField(@NotNull ReportableVariant variant) {
        String consequence = !variant.canonicalHgvsProteinImpact().isEmpty()
                ? AminoAcidFunctions.forceSingleLetterProteinAnnotation(variant.canonicalHgvsProteinImpact())
                : variant.canonicalHgvsCodingImpact();
        return variant.gene() + " " + consequence;
    }

    @NotNull
    private static String hotspotField(@NotNull ReportableVariant variant) {
        switch (variant.hotspot()) {
            case HOTSPOT:
                return "Yes";
            case NEAR_HOTSPOT:
                return "Near";
            default:
                return "No";
        }
    }

    @NotNull
    private static String driverLikelihoodField(@NotNull ReportableVariant variant) {
        if (Double.isNaN(variant.driverLikelihood())) {
            return Strings.EMPTY;
        } else {
            return PERCENTAGE_FORMAT.format(variant.driverLikelihood() * 100);
        }
    }
}
