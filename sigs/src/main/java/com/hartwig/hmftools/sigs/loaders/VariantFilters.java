package com.hartwig.hmftools.sigs.loaders;

import com.hartwig.hmftools.common.variant.SmallVariant;

import org.apache.commons.cli.CommandLine;

public class VariantFilters
{
    public final Double PloidyMin;
    public final Double PloidyMax;
    public final Double SubclonalLikelihoodMin;
    public final Double SubclonalLikelihoodMax;

    public static final String SUBCLONAL_MIN = "subclonal_min";
    public static final String SUBCLONAL_MAX = "subclonal_max";
    public static final String PLOIDY_MAX = "ploidy_max";
    public static final String PLOIDY_MIN = "ploidy_min";

    public VariantFilters(final CommandLine cmd)
    {
        SubclonalLikelihoodMin = initialiseDoubleValue(cmd, SUBCLONAL_MIN);
        SubclonalLikelihoodMax = initialiseDoubleValue(cmd, SUBCLONAL_MAX);

        PloidyMin = initialiseDoubleValue(cmd, PLOIDY_MIN);
        PloidyMax = initialiseDoubleValue(cmd, PLOIDY_MAX);
    }

    private static Double initialiseDoubleValue(final CommandLine cmd, final String config)
    {
        if(cmd.hasOption(config))
            return Double.parseDouble(cmd.getOptionValue(config, "0"));
        else
            return null;
    }

    public boolean passesFilters(final SmallVariant variant)
    {
        if(SubclonalLikelihoodMin != null && variant.subclonalLikelihood() < SubclonalLikelihoodMin)
            return false;

        if(SubclonalLikelihoodMax != null && variant.subclonalLikelihood() > SubclonalLikelihoodMax)
            return false;

        if(PloidyMin != null && variant.variantCopyNumber() < PloidyMin)
            return false;

        if(PloidyMax != null && variant.variantCopyNumber() > PloidyMax)
            return false;

        return true;
    }

}
