package com.hartwig.hmftools.sage.filter;

import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_VAF;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_GERMLINE_DEPTH;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MIN_TUMOR_VAF;
import static com.hartwig.hmftools.sage.filter.SoftFilter.MAX_GERMLINE_RELATIVE_QUAL;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.sage.common.VariantTier;

public class SoftFilterConfig
{
    public final String Name;
    public final double QualPScore;
    public final double MapQualFactor;
    public final double MinTumorVaf;
    public final int MinGermlineCoverage;
    public final int MinGermlineCoverageAllosome;
    public final double MaxGermlineVaf;
    public final double BaseMaxGermlineRelativeQual;


    private static final String MAP_QUAL_FACTOR = "map_qual_factor";

    public SoftFilterConfig(final ConfigBuilder configBuilder, final String prefix, final SoftFilterConfig defaultValue)
    {
        Name = defaultValue.Name;
        QualPScore = configBuilder.getDecimal(prefix + "_" + MIN_TUMOR_QUAL.configName());
        MapQualFactor = configBuilder.getDecimal(prefix + "_" + MAP_QUAL_FACTOR);
        MinTumorVaf = configBuilder.getDecimal(prefix + "_" + MIN_TUMOR_VAF.configName());
        MinGermlineCoverage = configBuilder.getInteger(prefix + "_" + MIN_GERMLINE_DEPTH.configName());
        MinGermlineCoverageAllosome = defaultValue.MinGermlineCoverageAllosome;
        MaxGermlineVaf = configBuilder.getDecimal(prefix + "_" + MAX_GERMLINE_VAF.configName());
        BaseMaxGermlineRelativeQual = configBuilder.getDecimal(prefix + "_" + MAX_GERMLINE_RELATIVE_QUAL.configName());
    }

    public SoftFilterConfig(
            final String name, final double qualPScore, final int mapQualFactor, final double minTumorVaf, final int minGermlineCoverage,
            final int minGermlineCoverageAllosome, final double maxGermlineVaf, final double maxGermlineRelativeQual)
    {
        Name = name;
        QualPScore = qualPScore;
        MapQualFactor = mapQualFactor;
        MinTumorVaf = minTumorVaf;
        MinGermlineCoverage = minGermlineCoverage;
        MinGermlineCoverageAllosome = minGermlineCoverageAllosome;
        MaxGermlineVaf = maxGermlineVaf;
        BaseMaxGermlineRelativeQual = maxGermlineRelativeQual;
    }

    public static SoftFilterConfig getTieredSoftFilterConfig(final VariantTier tier, final FilterConfig filterConfig)
    {
        switch(tier)
        {
            case HOTSPOT:
                return filterConfig.SoftHotspotFilter;
            case PANEL:
                return filterConfig.SoftPanelFilter;
            case HIGH_CONFIDENCE:
                return filterConfig.SoftHighConfidenceFilter;
            default:
                return filterConfig.SoftLowConfidenceFilter;
        }
    }

    public static void registerConfig(final ConfigBuilder configBuilder, final SoftFilterConfig defaultConfig)
    {
        String prefix = defaultConfig.Name;

        configBuilder.addDecimal(
                prefix + "_" + MIN_TUMOR_QUAL.configName(), "Minimum " + prefix + " tumor quality P-score", defaultConfig.QualPScore);

        configBuilder.addDecimal(
                prefix + "_" + MAP_QUAL_FACTOR, "Minimum " + prefix + " map qual factor",
                defaultConfig.MapQualFactor);

        configBuilder.addDecimal(
                prefix + "_" + MIN_TUMOR_VAF.configName(), "Minimum " + prefix + " tumor VAF",defaultConfig.MinTumorVaf);

        configBuilder.addInteger(
                prefix + "_" + MIN_GERMLINE_DEPTH.configName(),
                "Minimum " + prefix + " germline depth", defaultConfig.MinGermlineCoverage);

        configBuilder.addDecimal(
                prefix + "_" + MAX_GERMLINE_VAF.configName(), "Maximum " + prefix + " germline VAF", defaultConfig.MaxGermlineVaf);

        configBuilder.addDecimal(
                prefix + "_" + MAX_GERMLINE_RELATIVE_QUAL.configName(), "Base maximum " + prefix + " relative germline qual", defaultConfig.BaseMaxGermlineRelativeQual);
    }
}
