package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.hartwig.hmftools.orange.report.util.VariantUtil;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class GermlineVariantTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");

    private GermlineVariantTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, @NotNull List<ReportableVariant> variants) {
        if (variants.isEmpty()) {
            return TableUtil.createEmptyTable(title);
        }

        Table table = TableUtil.createReportContentTable(new float[] { 3, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Variant"), TableUtil.createHeaderCell("CN"), TableUtil.createHeaderCell("MACN"),
                        TableUtil.createHeaderCell("VCN"), TableUtil.createHeaderCell("RNA VAF"), TableUtil.createHeaderCell("Biallelic"),
                        TableUtil.createHeaderCell("Hotspot"), TableUtil.createHeaderCell("Genotype") });

        for (ReportableVariant variant : VariantUtil.sort(variants)) {
            table.addCell(TableUtil.createContentCell(VariantUtil.variantField(variant)));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(variant.totalCopyNumber())));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(variant.minorAlleleCopyNumber())));
            table.addCell(TableUtil.createContentCell(SINGLE_DIGIT.format(variant.alleleCopyNumber())));
            table.addCell(TableUtil.createContentCell("NA"));
            table.addCell(TableUtil.createContentCell(variant.biallelic() ? "Yes" : "No"));
            table.addCell(TableUtil.createContentCell(VariantUtil.hotspotField(variant)));
            table.addCell(TableUtil.createContentCell(variant.genotypeStatus().simplifiedDisplay()));
        }

        return TableUtil.createWrappingReportTable(table, title);
    }
}
