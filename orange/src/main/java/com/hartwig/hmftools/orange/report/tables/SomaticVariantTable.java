package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.interpretation.Variants;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class SomaticVariantTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<VariantEntry> variants,
            @NotNull ReportResources reportResources)
    {
        if(variants.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 3, 1, 1, 1, 1, 1, 1, 1, 1, 2 },
                new Cell[] { cells.createHeader("Variant"), cells.createHeader("VCN"), cells.createHeader("CN"), cells.createHeader("MACN"),
                        cells.createHeader("Biallelic"), cells.createHeader("Hotspot"), cells.createHeader("DL"), cells.createHeader("CL"),
                        cells.createHeader("Phase ID"), cells.createHeader("RNA Depth") });

        for(VariantEntry variant : Variants.sort(variants))
        {
            table.addCell(cells.createContent(Variants.variantField(variant)));
            table.addCell(cells.createContent(formatSingleDigitDecimal(variant.variantCopyNumber())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(variant.totalCopyNumber())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(variant.minorAlleleCopyNumber())));
            table.addCell(cells.createContent(Variants.biallelicLikelihoodField(variant)));
            table.addCell(cells.createContent(Variants.hotspotField(variant)));
            table.addCell(cells.createContent(Variants.driverLikelihoodField(variant)));
            table.addCell(cells.createContent(Variants.clonalLikelihoodField(variant)));
            table.addCell(cells.createContent(Variants.phaseSetField(variant)));
            table.addCell(cells.createContent(Variants.rnaDepthField(variant)));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }
}
