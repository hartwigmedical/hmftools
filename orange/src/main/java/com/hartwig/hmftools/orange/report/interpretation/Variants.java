package com.hartwig.hmftools.orange.report.interpretation;

import static java.lang.String.format;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.common.AllelicDepth;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.VariantEntry;

import org.apache.logging.log4j.util.Strings;

public final class Variants
{
    public static final String COL_VCN = "VCN";
    public static final String COL_MACN = "MACN";
    public static final String COL_CL = "Clonal";
    public static final String COL_DRIVER = "Driver";
    public static final String COL_AF = "AF";
    public static final String COL_SL = "Somatic";
    public static final String COL_VARIANT = "Variant";
    public static final String COL_BIALLELIC = "Biallelic";
    public static final String COL_HOTSPOT = "Hotspot";

    public static List<VariantEntry> sort(final List<VariantEntry> variants)
    {
        return variants.stream().sorted((variant1, variant2) ->
        {
            int driverCompare = Double.compare(variant2.driverLikelihood(), variant1.driverLikelihood());
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

    public static String rnaInfoField(final VariantEntry variant)
    {
        AllelicDepth rnaDepth = variant.rnaDepth();

        if(rnaDepth == null)
        {
            return ReportResources.NOT_AVAILABLE;
        }

        return format("%d/%d", rnaDepth.alleleReadCount(), rnaDepth.totalReadCount());
    }
}
