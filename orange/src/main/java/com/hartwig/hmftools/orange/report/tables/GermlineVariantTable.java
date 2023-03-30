package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;

import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.interpretation.Variants;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class GermlineVariantTable {

    private GermlineVariantTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<VariantEntry> variants) {
        if (variants.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 3, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Variant"), Cells.createHeader("VCN"), Cells.createHeader("CN"), Cells.createHeader("MACN"),
                        Cells.createHeader("RNA Depth"), Cells.createHeader("Biallelic"), Cells.createHeader("Hotspot"),
                        Cells.createHeader("Genotype") });

        for (VariantEntry variant : Variants.sort(variants)) {
            table.addCell(Cells.createContent(Variants.variantField(variant)));
            table.addCell(Cells.createContent(formatSingleDigitDecimal(variant.variantCopyNumber())));
            table.addCell(Cells.createContent(formatSingleDigitDecimal(variant.totalCopyNumber())));
            table.addCell(Cells.createContent(formatSingleDigitDecimal(variant.minorAlleleCopyNumber())));
            table.addCell(Cells.createContent(Variants.rnaDepthField(variant)));
            table.addCell(Cells.createContent(variant.biallelic() ? "Yes" : "No"));
            table.addCell(Cells.createContent(Variants.hotspotField(variant)));
            table.addCell(Cells.createContent(simplifiedDisplay(variant.genotypeStatus())));
        }

        return Tables.createWrapping(table, title);
    }

    private static String simplifiedDisplay(PurpleGenotypeStatus genotypeStatus) {
        switch (genotypeStatus) {
            case HOM_REF:
            case HOM_ALT:
                return "HOM";
            case HET: return "HET";
            case UNKNOWN: return "UNKNOWN";
        }
        throw new IllegalStateException();
    }
}
