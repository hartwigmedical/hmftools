package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;
import static com.hartwig.hmftools.orange.report.ReportResources.formatTwoDigitDecimal;

import java.text.DecimalFormat;
import java.util.List;

import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.itextpdf.layout.element.Cell;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.Nullable;

public final class TableCommon
{
    public static final String COL_GENE = "Gene";
    public static final String COL_RNA = "RNA";
    public static final String COL_TPM = "TPM";

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
