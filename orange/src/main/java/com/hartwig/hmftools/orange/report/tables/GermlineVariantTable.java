package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTwoDigitDecimal;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_AF;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_BIALLELIC;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_CN;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_DP;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_GENOTYPE;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_HOTSPOT;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_MACN;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_VARIANT;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_VCN;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.addEntry;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.cellArray;
import static com.hartwig.hmftools.orange.report.tables.SomaticVariantTable.floatArray;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.interpretation.Variants;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

public final class GermlineVariantTable
{
    public static Table build(final String title, float width, final List<VariantEntry> variants, final ReportResources reportResources)
    {
        if(variants.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 3, COL_VARIANT);
        addEntry(cells, widths, cellEntries, 1, COL_AF);
        addEntry(cells, widths, cellEntries, 1, COL_DP);
        addEntry(cells, widths, cellEntries, 1, COL_VCN);
        addEntry(cells, widths, cellEntries, 1, COL_CN);
        addEntry(cells, widths, cellEntries, 1, COL_MACN);
        addEntry(cells, widths, cellEntries, 1, COL_BIALLELIC);
        addEntry(cells, widths, cellEntries, 1, COL_HOTSPOT);
        addEntry(cells, widths, cellEntries, 1, COL_GENOTYPE);

        float[] widthArray = floatArray(widths);
        Cell[] cellArray = cellArray(cellEntries);

        Table table = Tables.createContent(width, widthArray, cellArray);

        for(VariantEntry variant : Variants.sort(variants))
        {
            table.addCell(cells.createContent(Variants.variantField(variant)));
            table.addCell(cells.createContent(formatTwoDigitDecimal(variant.vaf())));
            table.addCell(cells.createContent(String.valueOf(variant.depth())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(variant.variantCopyNumber())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(variant.totalCopyNumber())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(variant.minorAlleleCopyNumber())));
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
