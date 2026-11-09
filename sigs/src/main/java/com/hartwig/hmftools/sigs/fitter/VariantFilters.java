package com.hartwig.hmftools.sigs.fitter;

import static com.hartwig.hmftools.common.variant.PaveVcfTags.GNOMAD_FREQ;
import static com.hartwig.hmftools.sigs.fitter.SigsConfig.TUMOR_ONLY_AF_THRESHOLD;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.common.variant.SomaticLikelihood;

public class VariantFilters
{
    public final Double PloidyMin;
    public final Double PloidyMax;
    public final Double SubclonalLikelihoodMin;
    public final Double SubclonalLikelihoodMax;
    public final boolean TumorOnly;

    public static final String SUBCLONAL_MIN = "subclonal_min";
    public static final String SUBCLONAL_MAX = "subclonal_max";
    public static final String PLOIDY_MAX = "ploidy_max";
    public static final String PLOIDY_MIN = "ploidy_min";
    public static final String TUMOR_ONLY = "tumor_only";

    public VariantFilters(final ConfigBuilder configBuilder)
    {
        SubclonalLikelihoodMin = initialiseDoubleValue(configBuilder, SUBCLONAL_MIN);
        SubclonalLikelihoodMax = initialiseDoubleValue(configBuilder, SUBCLONAL_MAX);

        PloidyMin = initialiseDoubleValue(configBuilder, PLOIDY_MIN);
        PloidyMax = initialiseDoubleValue(configBuilder, PLOIDY_MAX);

        TumorOnly = configBuilder.hasFlag(TUMOR_ONLY);
    }

    private static Double initialiseDoubleValue(final ConfigBuilder configBuilder, final String configStr)
    {
        if(configBuilder.hasValue(configStr))
            return configBuilder.getDecimal(configStr);
        else
            return null;
    }

    public boolean passesFilters(final SigVariant variant)
    {
        SmallVariant smallVariant = variant.Variant;

        if(SubclonalLikelihoodMin != null && smallVariant.subclonalLikelihood() < SubclonalLikelihoodMin)
            return false;

        if(SubclonalLikelihoodMax != null && smallVariant.subclonalLikelihood() > SubclonalLikelihoodMax)
            return false;

        if(PloidyMin != null && smallVariant.variantCopyNumber() < PloidyMin)
            return false;

        if(PloidyMax != null && smallVariant.variantCopyNumber() > PloidyMax)
            return false;

        if(TumorOnly)
        {
            if(smallVariant.somaticLikelihood() == SomaticLikelihood.LOW)
                return false;

            if(smallVariant.isHotspot())
                return false;

            // non-zero GND_FREQ
            double gnomadFreq = variant.VariantContext.getAttributeAsDouble(GNOMAD_FREQ, -1);
            if(gnomadFreq > 0)
                return false;

            // AF < 5%
            double af = variant.Variant.allelicDepth().alleleFrequency();;
            if(af < TUMOR_ONLY_AF_THRESHOLD)
                return false;
        }

        return true;
    }

    public static void registerConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addDecimal(SUBCLONAL_MIN,"Subclonal minimum", 0);
        configBuilder.addDecimal(SUBCLONAL_MAX,"Subclonal minimum", 1.0);
        configBuilder.addDecimal(PLOIDY_MIN,"Ploidy minimum", 0);
        configBuilder.addDecimal(PLOIDY_MAX,"Ploidy minimum", 1.0);
        configBuilder.addFlag(TUMOR_ONLY, "Apply tumor-only filters");
    }

    @VisibleForTesting
    public VariantFilters(
            final Double ploidyMin, final Double ploidyMax, final Double subclonalLikelihoodMin,
            final Double subclonalLikelihoodMax, final boolean tumorOnly)
    {
        PloidyMin = ploidyMin;
        PloidyMax = ploidyMax;
        SubclonalLikelihoodMin = subclonalLikelihoodMin;
        SubclonalLikelihoodMax = subclonalLikelihoodMax;
        TumorOnly = tumorOnly;
    }
}
