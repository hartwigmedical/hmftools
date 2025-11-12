package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;

import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.datamodel.finding.SmallVariant;
import com.hartwig.hmftools.datamodel.finding.Variants;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class GermlineVariantTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<SmallVariant> variants,
            @NotNull ReportResources reportResources)
    {
        if(variants.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 3, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { cells.createHeader("Variant"), cells.createHeader("VCN"), cells.createHeader("CN"), cells.createHeader("MACN"),
                        cells.createHeader("RNA Depth"), cells.createHeader("Biallelic"), cells.createHeader("Hotspot"),
                        cells.createHeader("Genotype") });

        for(SmallVariant variant : Variants.sort(variants))
        {
            table.addCell(cells.createContent(Variants.variantField(variant, true)));
            table.addCell(cells.createContent(formatSingleDigitDecimal(variant.variantCopyNumber())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(variant.totalCopyNumber())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(variant.minorAlleleCopyNumber())));
            table.addCell(cells.createContent(Variants.rnaDepthField(variant, ReportResources.NOT_AVAILABLE)));
            table.addCell(cells.createContent(variant.biallelic() ? "Yes" : "No"));
            table.addCell(cells.createContent(Variants.hotspotField(variant)));
            table.addCell(cells.createContent(simplifiedDisplay(variant.genotypeStatus())));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static String simplifiedDisplay(PurpleGenotypeStatus genotypeStatus)
    {
        switch(genotypeStatus)
        {
            case HOM_REF:
            case HOM_ALT:
                return "HOM";
            case HET:
                return "HET";
            case UNKNOWN:
                return "UNKNOWN";
        }
        throw new IllegalStateException();
    }
}
