package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.hartwig.hmftools.orange.report.util.VariantUtil;
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
    public static Table build(@NotNull String title, float width, @NotNull List<ReportableVariant> driverVariants) {
        if (driverVariants.isEmpty()) {
            return TableUtil.createEmptyTable(title, width);
        }

        Table table = TableUtil.createReportContentTable(width,
                new float[] { 3, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Variant"), TableUtil.createHeaderCell("CN"), TableUtil.createHeaderCell("MACN"),
                        TableUtil.createHeaderCell("VCN"), TableUtil.createHeaderCell("RNA VAF"), TableUtil.createHeaderCell("Biallelic"),
                        TableUtil.createHeaderCell("Hotspot"), TableUtil.createHeaderCell("DL"), TableUtil.createHeaderCell("CL"),
                        TableUtil.createHeaderCell("Phase") });

        for (ReportableVariant variant : VariantUtil.sort(driverVariants)) {
            table.addCell(TableUtil.createContentCell(VariantUtil.variantField(variant)));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(variant.totalCopyNumber())));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(variant.minorAlleleCopyNumber())));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(variant.alleleCopyNumber())));
            table.addCell(TableUtil.createContentCell(ReportResources.NOT_AVAILABLE));
            table.addCell(TableUtil.createContentCell(variant.biallelic() ? "Yes" : "No"));
            table.addCell(TableUtil.createContentCell(VariantUtil.hotspotField(variant)));
            table.addCell(TableUtil.createContentCell(driverLikelihoodField(variant)));
            table.addCell(TableUtil.createContentCell(PERCENTAGE_FORMAT.format(variant.clonalLikelihood() * 100)));
            table.addCell(TableUtil.createContentCell(
                    variant.localPhaseSet() != null ? String.valueOf(variant.localPhaseSet()) : Strings.EMPTY));
        }

        return TableUtil.createWrappingReportTable(table, title);
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
