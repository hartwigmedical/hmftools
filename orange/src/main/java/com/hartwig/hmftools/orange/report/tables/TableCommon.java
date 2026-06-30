package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

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

    public static final String COL_SUPPORT = "Support";
    public static final String COL_SUPPORT_START = "Support Start";
    public static final String COL_SUPPORT_END = "Support End";

    public static final String COL_RNA = "RNA";
    public static final String COL_RNA_FRAGS = "RNA Frags";
    public static final String COL_TPM = "TPM";

    public static final String COL_JUNC_START = "Junc Start";
    public static final String COL_JUNC_END = "Junc End";

    public static final String VALUE_HET = "HET";
    public static final String VALUE_HOM = "HOM";

    public static final String VALUE_YES = "Yes";
    public static final String VALUE_NO = "No";

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

    private static String formatDecimal(double num, final String fmt)
    {
        return new DecimalFormat(fmt, DecimalFormatSymbols.getInstance(Locale.ENGLISH)).format(num);
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
            return location;
        }
    }
}
