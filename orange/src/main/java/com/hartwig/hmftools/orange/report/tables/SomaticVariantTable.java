package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.interpretation.Variants;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class SomaticVariantTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");

    private SomaticVariantTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<VariantEntry> variants,
            @NotNull ReportResources reportResources) {
        if (variants.isEmpty()) {
            return reportResources.tables().createEmpty(title, width);
        }

        Cells cells = reportResources.cells();
        Table table = Tables.createContent(width,
                new float[] { 3, 1, 1, 1, 1, 1, 1, 1, 1, 2 },
                new Cell[] { cells.createHeader("Variant"), cells.createHeader("VCN"), cells.createHeader("CN"), cells.createHeader("MACN"),
                        cells.createHeader("Biallelic"), cells.createHeader("Hotspot"), cells.createHeader("DL"), cells.createHeader("CL"),
                        cells.createHeader("Phase ID"), cells.createHeader("RNA Depth") });

        for (VariantEntry variant : Variants.sort(variants)) {
            table.addCell(cells.createContent(Variants.variantField(variant)));
            table.addCell(cells.createContent(SINGLE_DIGIT.format(variant.variantCopyNumber())));
            table.addCell(cells.createContent(SINGLE_DIGIT.format(variant.totalCopyNumber())));
            table.addCell(cells.createContent(SINGLE_DIGIT.format(variant.minorAlleleCopyNumber())));
            table.addCell(cells.createContent(variant.biallelic() ? "Yes" : "No"));
            table.addCell(cells.createContent(Variants.hotspotField(variant)));
            table.addCell(cells.createContent(Variants.driverLikelihoodField(variant)));
            table.addCell(cells.createContent(Variants.clonalLikelihoodField(variant)));
            table.addCell(cells.createContent(Variants.phaseSetField(variant)));
            table.addCell(cells.createContent(Variants.rnaDepthField(variant)));
        }

        return reportResources.tables().createWrapping(table, title);
    }
}
