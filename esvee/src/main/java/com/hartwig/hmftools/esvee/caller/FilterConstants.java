package com.hartwig.hmftools.esvee.caller;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V37;
import static com.hartwig.hmftools.common.region.ExcludedRegions.getPolyGRegions;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TARGET_REGIONS_BED;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.DECIMAL;
import static com.hartwig.hmftools.common.utils.config.ConfigItemType.INTEGER;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_VARIANT_LENGTH;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class FilterConstants
{
    public final int MinQual;

    public static final int DEFAULT_MIN_TUMOR_QUAL = 30;

    public final int MinLength;

    public final int MinSupportJunction;
    public final int MinSupportSgl;
    public final int MinSupportHotspot;

    public static final int DEFAULT_MIN_SUPPORT_JUNCTION = 4;
    public static final int DEFAULT_MIN_SUPPORT_SGL = 6;
    public static final int DEFAULT_MIN_SUPPORT_HOTSPOT = 2;

    public final double MinAfJunction;
    public final double MinAfSgl;
    public final double MinAfHotspot;

    public static final double DEFAULT_MIN_AF_JUNCTION = 0.01;
    public static final double DEFAULT_MIN_AF_SGL = 0.05;
    public static final double DEFAULT_MIN_AF_HOTSPOT = 0.001;

    public final int MinAvgFragFactor;
    public static final int DEFAULT_MIN_AVG_FRAG_FACTOR = 3;

    public final boolean FilterSGLs;
    private static final String FILTER_SGLS = "filter_sgls";

    public final List<ChrBaseRegion> PolyGcRegions;
    public final ChrBaseRegion LowQualRegion;

    public final int PonDistance;
    public static final int DEFAULT_PON_DISTANCE = 3;
    public static final String PON_DISTANCE = "pon_distance";

    public static final int SHORT_CALLING_SIZE = 1000;
    public static final double SGL_MIN_STRAND_BIAS = 0.05;
    public static final double SGL_MAX_STRAND_BIAS = 0.95;
    public static final double MAX_STRAND_BIAS = 0.95;

    public static final ChrBaseRegion PMS2_V37 = new ChrBaseRegion("7", 6002870, 6058756); // has 10K buffer
    public static final ChrBaseRegion PMS2_V38 = new ChrBaseRegion("chr7", 5960925, 6019106);

    public static final double GERMLINE_AF_THRESHOLD = 0.1;

    public static FilterConstants from(final ConfigBuilder configBuilder)
    {
        boolean targetedMode = configBuilder.hasValue(TARGET_REGIONS_BED);
        RefGenomeVersion refGenVersion = RefGenomeVersion.from(configBuilder);

        // use targeted panel defaults where applicable
        boolean filterSgls = targetedMode || configBuilder.hasFlag(FILTER_SGLS);

        return new FilterConstants(
                DEFAULT_MIN_TUMOR_QUAL, MIN_VARIANT_LENGTH,
                DEFAULT_MIN_SUPPORT_JUNCTION, DEFAULT_MIN_SUPPORT_SGL, DEFAULT_MIN_SUPPORT_HOTSPOT,
                DEFAULT_MIN_AF_JUNCTION, DEFAULT_MIN_AF_SGL, DEFAULT_MIN_AF_HOTSPOT,
                DEFAULT_MIN_AVG_FRAG_FACTOR, filterSgls,
                getPolyGRegions(refGenVersion), refGenVersion == V37 ? PMS2_V37 : PMS2_V38,
                configBuilder.getInteger(PON_DISTANCE));
    }

    public FilterConstants(
            final int minTumorQual, final int minLength, final int minSupportJunction, final int minSupportSgl, final int minSupportHotspot,
            final double minAfJunction, final double minAfSgl, final double minAfHotspot, final int minAvgFragFactor,
            final boolean filterSGLs, final List<ChrBaseRegion> polyGcRegions, final ChrBaseRegion lowQualRegion, final int ponDistance)
    {
        MinQual = minTumorQual;
        MinLength = minLength;
        MinSupportJunction = minSupportJunction;
        MinSupportSgl = minSupportSgl;
        MinSupportHotspot = minSupportHotspot;
        MinAfJunction = minAfJunction;
        MinAfSgl = minAfSgl;
        MinAfHotspot = minAfHotspot;
        MinAvgFragFactor = minAvgFragFactor;
        FilterSGLs = filterSGLs;
        PolyGcRegions = polyGcRegions;
        LowQualRegion = lowQualRegion;
        PonDistance = ponDistance;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        configBuilder.addInteger(PON_DISTANCE, "PON permitted margin", DEFAULT_PON_DISTANCE);
        configBuilder.addFlag(FILTER_SGLS, "Filter SGLs from VCF, intended for tumor-only mode, default=true in target panel");

        /*
        addTargetedInteger(
                configBuilder, HARD_MIN_TUMOR_QUAL_CFG, "Hard min tumor qual",
                DEFAULT_HARD_MIN_TUMOR_QUAL, TARGETED_DEFAULT_HARD_MIN_TUMOR_QUAL);

        addTargetedInteger(
                configBuilder, MIN_QUAL_BREAK_END_CFG, "Min qual break end",
                DEFAULT_MIN_QUAL_BREAK_END, TARGETED_DEFAULT_MIN_QUAL_BREAK_END);

        addTargetedInteger(
                configBuilder, MIN_QUAL_BREAK_POINT_CFG, "Min qual break point",
                DEFAULT_MIN_QUAL_BREAK_POINT, TARGETED_DEFAULT_MIN_QUAL_BREAK_POINT);

        addTargetedInteger(configBuilder, QUAL_PER_AD, "Qual per AD limit", 0, TARGETED_DEFAULT_QUAL_PER_AD);

        addTargetedDecimal(configBuilder, MODIFIED_AF, "Modified AF limit", 0, TARGETED_DEFAULT_MODIFIED_AF);

        addTargetedDecimal(
                configBuilder, MODIFIED_AF_HOTSPOT, "Modified AF limit for hotspots", 0, TARGETED_DEFAULT_MODIFIED_AF_HOTSPOT);
        */

    }

    private static void addTargetedDecimal(
            final ConfigBuilder configBuilder, final String name, final String desc, double defaultValue, double targetedDefaultValue)
    {
        configBuilder.addConfigItem(
                DECIMAL, name, false,
                format("%s, default=%.3g targeted default=%.3g", desc, defaultValue, targetedDefaultValue), String.valueOf(defaultValue));
    }

    private static void addTargetedInteger(
            final ConfigBuilder configBuilder, final String name, final String desc, int defaultValue, int targetedDefaultValue)
    {
        configBuilder.addConfigItem(
                INTEGER, name, false,
                format("%s, default=%d targeted default=%d", desc, defaultValue, targetedDefaultValue), String.valueOf(defaultValue));
    }
}
