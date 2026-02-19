package com.hartwig.hmftools.orange.report.interpretation;

import static java.lang.String.format;

import java.text.DecimalFormat;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class Variants
{
    private static final DecimalFormat PERCENTAGE_FORMAT = new DecimalFormat("#'%'");

    public static final String COL_VCN = "VCN";
    public static final String COL_CN = "CN";
    public static final String COL_MACN = "MACN";
    public static final String COL_CL = "CL";
    public static final String COL_DL = "DL";
    public static final String COL_AF = "AF";
    public static final String COL_DP = "DEPTH";
    public static final String COL_SL = "SL";
    public static final String COL_VARIANT = "Variant";
    public static final String COL_GENOTYPE = "Genotype";
    public static final String COL_BIALLELIC = "Biallelic";
    public static final String COL_HOTSPOT = "Hotspot";
    public static final String COL_RNA = "RNA";

    public static List<VariantEntry> sort(final List<VariantEntry> variants)
    {
        return variants.stream().sorted((variant1, variant2) ->
        {
            double driverLikelihood1 = variant1.driverLikelihood() != null ? variant1.driverLikelihood() : -1;
            double driverLikelihood2 = variant2.driverLikelihood() != null ? variant2.driverLikelihood() : -1;

            int driverCompare = Double.compare(driverLikelihood2, driverLikelihood1);
            if(driverCompare != 0)
            {
                return driverCompare;
            }

            int geneCompare = variant1.gene().compareTo(variant2.gene());
            if(geneCompare != 0)
            {
                return geneCompare;
            }

            if(variant1.affectedCodon() == null && variant2.affectedCodon() == null)
            {
                return 0;
            }
            else if(variant1.affectedCodon() == null)
            {
                return 1;
            }
            else if(variant2.affectedCodon() == null)
            {
                return -1;
            }
            else
            {
                return Integer.compare(variant1.affectedCodon(), variant2.affectedCodon());
            }
        }).collect(Collectors.toList());
    }

    public static String variantField(final VariantEntry variant)
    {
        String addon = Strings.EMPTY;

        if(!variant.isCanonical())
        {
            addon = " (alt)";
        }

        return variant.gene() + addon + " " + variant.impact();
    }

    public static String hotspotField(final VariantEntry variant)
    {
        switch(variant.hotspot())
        {
            case HOTSPOT:
                return "Yes";
            case NEAR_HOTSPOT:
                return "Near";
            default:
                return "No";
        }
    }

    public static String biallelicLikelihoodField(final VariantEntry variant)
    {
        return PERCENTAGE_FORMAT.format(variant.biallelicProbability() * 100);
    }

    public static String driverLikelihoodField(final VariantEntry variant)
    {
        return variant.driverLikelihood() != null ? PERCENTAGE_FORMAT.format(variant.driverLikelihood() * 100) : Strings.EMPTY;
    }

    public static String clonalLikelihoodField(final VariantEntry variant)
    {
        return PERCENTAGE_FORMAT.format(100 * variant.clonalLikelihood());
    }

    public static String rnaInfoField(final VariantEntry variant)
    {
        PurpleAllelicDepth rnaDepth = variant.rnaDepth();

        if(rnaDepth == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }

        return format("%d/%d", rnaDepth.alleleReadCount(), rnaDepth.totalReadCount());
    }
}
