package com.hartwig.hmftools.orange.report.tables;

import java.util.List;

import com.hartwig.hmftools.orange.report.util.Cells;
import com.itextpdf.layout.element.Cell;

public final class TableCommon
{
    public static final String COL_GENE = "Gene";
    public static final String COL_CHR = "Chromosome";
    public static final String COL_CN = "CN";
    public static final String COL_JCN = "JCN";
    public static final String COL_REL_CN = "Rel CN";
    public static final String COL_RNA = "RNA";
    public static final String COL_TPM = "TPM";
    public static final String COL_LOCATION = "Location";
    public static final String COL_RANGE = "Range";
    public static final String COL_TYPE = "Type";
    public static final String COL_DRIVER = "Driver";
    public static final String COL_ZYGOSITY = "Zygosity";

    public static final String VALUE_HET = "HET";
    public static final String VALUE_HOM = "HOM";

    protected static void addEntry(final Cells cells, final List<Integer> widths, final List<Cell> cellEntries, int width, final String column)
    {
        cellEntries.add(cells.createHeader(column));
        widths.add(width);
    }

    protected static void addEntry(final Cells cells, final List<Float> widths, final List<Cell> cellEntries, double width, final String column)
    {
        cellEntries.add(cells.createHeader(column));
        widths.add((float)width);
    }

    protected static float[] intToFloatArray(final List<Integer> widths)
    {
        float[] widthArray = new float[widths.size()];

        for(int i = 0; i < widthArray.length; ++i)
        {
            widthArray[i] = widths.get(i);
        }

        return widthArray;
    }

    protected static float[] floatArray(final List<Float> widths)
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
