package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTwoDigitDecimal;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_AF;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_BIALLELIC;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_CL;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_CN;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_DL;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_DP;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_HOTSPOT;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_MACN;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_RNA;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_SL;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_VARIANT;
import static com.hartwig.hmftools.orange.report.interpretation.Variants.COL_VCN;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.interpretation.Variants;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

public final class SomaticVariantTable
{
    public static Table build(
            final String title, float width, final List<VariantEntry> variants, final ReportResources reportResources,
            boolean tumorOnly, boolean hasRna)
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
        addEntry(cells, widths, cellEntries, 1, COL_DL);
        addEntry(cells, widths, cellEntries, 1, COL_CL);

        if(tumorOnly)
        {
            addEntry(cells, widths, cellEntries, 1, COL_SL);
        }

        if(hasRna)
        {
            addEntry(cells, widths, cellEntries, 1, COL_RNA);
        }

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
            table.addCell(cells.createContent(Variants.biallelicLikelihoodField(variant)));
            table.addCell(cells.createContent(Variants.hotspotField(variant)));
            table.addCell(cells.createContent(Variants.driverLikelihoodField(variant)));
            table.addCell(cells.createContent(Variants.clonalLikelihoodField(variant)));

            if(tumorOnly)
                table.addCell(cells.createContent(variant.somaticLikelihood()));

            if(hasRna)
                table.addCell(cells.createContent(Variants.rnaInfoField(variant)));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    protected static void addEntry(final Cells cells, final List<Integer> widths, final List<Cell> cellEntries, int width, final String column)
    {
        cellEntries.add(cells.createHeader(column));
        widths.add(width);
    }

    protected static float[] floatArray(final List<Integer> widths)
    {
        float[] widthArray = new float[widths.size()];

        for(int i = 0; i < widthArray.length; ++i)
        {
            widthArray[i] = widths.get(i);
        }

        return widthArray;
    }

    protected static Cell[] cellArray(final List<Cell> cellEntries)
    {
        Cell[] cellArray = new Cell[cellEntries.size()];

        for(int i = 0; i < cellArray.length; ++i)
        {
            cellArray[i] = cellEntries.get(i);
        }

        return cellArray;
    }
}

