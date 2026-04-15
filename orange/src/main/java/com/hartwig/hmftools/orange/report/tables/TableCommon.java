package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.List;
import java.util.Locale;

import java.io.IOException;

import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;

import org.apache.pdfbox.pdmodel.PDPage;

import be.quodlibet.boxable.BaseTable;
import be.quodlibet.boxable.Cell;
import be.quodlibet.boxable.Row;

public final class TableCommon
{
    public static final String COL_GENE = "Gene";
    public static final String COL_CHR = "Chromosome";
    public static final String COL_POSITION = "Position";
    public static final String COL_CN = "CN";
    public static final String COL_JCN = "JCN";
    public static final String COL_REL_CN = "Rel CN";
    public static final String COL_LOCATION = "Location";
    public static final String COL_JUNCTIONS = "Junctions";
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

    // variants
    public static final String COL_VARIANT = "Variant";
    public static final String COL_HGVS = "HGVS";
    public static final String COL_VCN = "VCN";
    public static final String COL_COPIES = "Copies";
    public static final String COL_MACN = "MACN";
    public static final String COL_CL = "Clonal";
    public static final String COL_AF = "AF";
    public static final String COL_SL = "Somatic";
    public static final String COL_BIALLELIC = "Biall";
    public static final String COL_HOTSPOT = "Hotspot";

    public static final String COL_SUPPORT = "Support"; // for fields of type fragments / depth
    public static final String COL_SUPPORT_START = "Support Start";
    public static final String COL_SUPPORT_END = "Support End";

    // RNA
    public static final String COL_RNA = "RNA";
    public static final String COL_RNA_FRAGS = "RNA Frags";
    public static final String COL_TPM = "TPM";

    public static final String COL_JUNC_START = "Junc Start";
    public static final String COL_JUNC_END = "Junc End";

    public static final String VALUE_HET = "HET";
    public static final String VALUE_HOM = "HOM";

    public static final String VALUE_YES = "Yes";
    public static final String VALUE_NO = "No";

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

    public static String formatTpmField(final double tpm)
    {
        return formatSingleDigitDecimal(tpm);
    }

    public static String formatPercentileField(final double percentile)
    {
        return formatTwoDigitDecimal(percentile);
    }

    public static String formatFoldChangeField(final double foldChange)
    {
        return foldChange > 1000 ? ">1000" : formatSingleDigitDecimal(foldChange);
    }

    public static String formatSupportField(final int fragments, final int depth)
    {
        return format("%d/%d", fragments, depth);
    }

    public static String zeroPrefixed(final String location)
    {
        int armStart = location.indexOf("q");
        if(armStart < 0)
        {
            armStart = location.indexOf("p");
        }

        String chromosome = armStart > 0 ? location.substring(0, armStart) : location;

        try
        {
            int chromosomeIndex = Integer.parseInt(chromosome);
            if(chromosomeIndex < 10)
            {
                return "0" + location;
            }
            else
            {
                return location;
            }
        }
        catch(NumberFormatException exception)
        {
            // No need to prefix Y/X chromosomes
            return location;
        }
    }

    // Column definition helpers - collect widths and header names
    public static void addEntry(final List<Integer> widths, final List<String> headers, int width, final String column)
    {
        widths.add(width);
        headers.add(column);
    }

    public static void addEntry(final List<Float> widths, final List<String> headers, double width, final String column)
    {
        widths.add((float) width);
        headers.add(column);
    }

    // Legacy overloads that still accept Cells parameter (ignored) for easier migration
    public static void addEntry(final Cells cells, final List<Integer> widths, final List<String> headers, int width, final String column)
    {
        addEntry(widths, headers, width, column);
    }

    public static void addEntry(final Cells cells, final List<Float> widths, final List<String> headers, double width, final String column)
    {
        addEntry(widths, headers, width, column);
    }

    public static float[] intToFloatArray(final List<Integer> widths)
    {
        float[] widthArray = new float[widths.size()];

        for(int i = 0; i < widthArray.length; ++i)
        {
            widthArray[i] = widths.get(i);
        }

        return widthArray;
    }

    public static float[] floatArray(final List<Float> widths)
    {
        float[] widthArray = new float[widths.size()];

        for(int i = 0; i < widthArray.length; ++i)
        {
            widthArray[i] = widths.get(i);
        }

        return widthArray;
    }

    public static String[] stringArray(final List<String> headers)
    {
        return headers.toArray(new String[0]);
    }

    public static BaseTable createStandardTable(
            final DocumentContext docCtx, final String title, float width,
            final float[] columnWidths, final String[] headerTexts,
            final ReportResources reportResources) throws IOException
    {
        return new Tables(reportResources).createWithTitle(docCtx, title, width, columnWidths, headerTexts);
    }

    public static BaseTable createEmptyTable(
            final DocumentContext docCtx, final String title, float width,
            final ReportResources reportResources) throws IOException
    {
        return new Tables(reportResources).createEmpty(docCtx, title, width);
    }

    public static float[] toPercentages(final float[] widths)
    {
        float total = 0;
        for(float w : widths)
        {
            total += w;
        }
        float[] pcts = new float[widths.length];
        for(int i = 0; i < widths.length; i++)
        {
            pcts[i] = (widths[i] / total) * 100f;
        }
        return pcts;
    }

    public static java.util.List<Cell<PDPage>> addDataRow(
            final BaseTable table, final Cells cells, final float[] pcts, final String[] values)
    {
        Row<PDPage> row = table.createRow(Tables.rowHeight());
        java.util.List<Cell<PDPage>> createdCells = new java.util.ArrayList<>();
        for(int i = 0; i < values.length; i++)
        {
            createdCells.add(cells.addContentCell(row, pcts[i], values[i]));
        }
        return createdCells;
    }
}
