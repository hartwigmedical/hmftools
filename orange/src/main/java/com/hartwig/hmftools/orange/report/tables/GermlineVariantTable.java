package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Variants;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class GermlineVariantTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");

    private GermlineVariantTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<ReportableVariant> variants) {
        if (variants.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 3, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Variant"), Cells.createHeader("VCN"), Cells.createHeader("CN"), Cells.createHeader("MACN"),
                        Cells.createHeader("RNA VAF"), Cells.createHeader("Biallelic"), Cells.createHeader("Hotspot"),
                        Cells.createHeader("Genotype") });

        for (ReportableVariant variant : Variants.sort(variants)) {
            table.addCell(Cells.createContent(Variants.variantField(variant)));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(variant.alleleCopyNumber())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(variant.totalCopyNumber())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(variant.minorAlleleCopyNumber())));
            table.addCell(Cells.createContent(ReportResources.NOT_AVAILABLE));
            table.addCell(Cells.createContent(variant.biallelic() ? "Yes" : "No"));
            table.addCell(Cells.createContent(Variants.hotspotField(variant)));
            table.addCell(Cells.createContent(variant.genotypeStatus().simplifiedDisplay()));
        }

        return Tables.createWrapping(table, title);
    }
}
