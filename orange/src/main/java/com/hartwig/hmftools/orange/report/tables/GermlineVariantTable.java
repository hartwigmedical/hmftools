package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.common.variant.ReportableVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.CellUtil;
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
    public static Table build(@NotNull String title, float width, @NotNull List<ReportableVariant> variants) {
        if (variants.isEmpty()) {
            return TableUtil.createEmpty(title, width);
        }

        Table table = TableUtil.createContent(width,
                new float[] { 3, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { CellUtil.createHeader("Variant"), CellUtil.createHeader("VCN"), CellUtil.createHeader("CN"),
                        CellUtil.createHeader("MACN"), CellUtil.createHeader("RNA VAF"), CellUtil.createHeader("Biallelic"),
                        CellUtil.createHeader("Hotspot"), CellUtil.createHeader("Genotype") });

        for (ReportableVariant variant : VariantUtil.sort(variants)) {
            table.addCell(CellUtil.createContent(VariantUtil.variantField(variant)));
            table.addCell(CellUtil.createContent(SINGLE_DIGIT.format(variant.alleleCopyNumber())));
            table.addCell(CellUtil.createContent(SINGLE_DIGIT.format(variant.totalCopyNumber())));
            table.addCell(CellUtil.createContent(SINGLE_DIGIT.format(variant.minorAlleleCopyNumber())));
            table.addCell(CellUtil.createContent(ReportResources.NOT_AVAILABLE));
            table.addCell(CellUtil.createContent(variant.biallelic() ? "Yes" : "No"));
            table.addCell(CellUtil.createContent(VariantUtil.hotspotField(variant)));
            table.addCell(CellUtil.createContent(variant.genotypeStatus().simplifiedDisplay()));
        }

        return TableUtil.createWrapping(table, title);
    }
}
