package com.hartwig.hmftools.sigs.fitter;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.SmallVariant;

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

    public VariantFilters(final ConfigBuilder configBuilder)
    {
        SubclonalLikelihoodMin = initialiseDoubleValue(configBuilder, SUBCLONAL_MIN);
        SubclonalLikelihoodMax = initialiseDoubleValue(configBuilder, SUBCLONAL_MAX);

        PloidyMin = initialiseDoubleValue(configBuilder, PLOIDY_MIN);
        PloidyMax = initialiseDoubleValue(configBuilder, PLOIDY_MAX);
    }

    private static Double initialiseDoubleValue(final ConfigBuilder configBuilder, final String configStr)
    {
        if(configBuilder.hasValue(configStr))
            return configBuilder.getDecimal(configStr);
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

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addDecimal(SUBCLONAL_MIN,"Subclonal minimum", 0);
        configBuilder.addDecimal(SUBCLONAL_MAX,"Subclonal minimum", 1.0);
        configBuilder.addDecimal(PLOIDY_MIN,"Ploidy minimum", 0);
        configBuilder.addDecimal(PLOIDY_MAX,"Ploidy minimum", 1.0);
    }
}
