package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;

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
    public static final String COL_SV_TYPE = "SV Type";
    public static final String COL_DRIVER = "Driver";
    public static final String COL_ZYGOSITY = "Zygosity";

    public static final String COL_COHOR_FREQ = "Cohort Freq";
    public static final String COL_FRAGS = "Fragments";
    public static final String COL_DEPTH = "Depth";
    public static final String COL_DP = COL_DEPTH;
    public static final String COL_DEPTH_START = "Depth Start";
    public static final String COL_DEPTH_END = "Depth End";
    public static final String COL_FUSION = "Fusion";

    public static final String COL_SUPPORT = "Support"; // for fields of type fragments / depth
    public static final String COL_SUPPORT_START = "Support Start";
    public static final String COL_SUPPORT_END = "Support End";

    public static final String COL_JUNC_START = "Junc Start";
    public static final String COL_JUNC_END = "Junc End";

    public static final String VALUE_HET = "HET";
    public static final String VALUE_HOM = "HOM";

    // common field formatting

    // formatting
    public static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");

    public static String formatSingleDigitDecimal(double num)
    {
        return formatDecimal(num, "0.0");
    }
    public static String formatTwoDigitDecimal(double num)
    {
        return formatDecimal(num, "0.00");
    }
    public static String formatPercentage(double num)
    {
        return formatPercentage(num, true);
    }

    public static String formatPercentage(double num, boolean multiplyBy100)
    {
        return formatDecimal(multiplyBy100 ? num * 100 : num, "0'%'");
    }

    public static String formatPercentageOneDecimal(double num)
    {
        return formatDecimal(num * 100, "0.0'%'");
    }

    private static String formatDecimal(double num, final String format)
    {
        // To make sure every decimal format uses a dot as separator rather than a comma.
        return new DecimalFormat(format, DecimalFormatSymbols.getInstance(Locale.ENGLISH)).format(num);
    }

    public static String formatPercentageField(final double value)
    {
        return PERCENTAGE_FORMAT.format(value * 100);
    }

    public static String formatTpmField(final double tpm) { return formatSingleDigitDecimal(tpm); }

    public static String formatPercentileField(final double percentile) { return formatTwoDigitDecimal(percentile); }

    public static String formatFoldChangeField(final double foldChange)
    {
        return foldChange > 1000 ? ">1000" : formatSingleDigitDecimal(foldChange);
    }

    public static String formatSupportField(final int fragments, final int depth)
    {
        return format("%d/%d", fragments, depth);
    }

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
