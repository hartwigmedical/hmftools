package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.sage.config.BaseQualityRecalibrationConfig.DEFAULT_BQR_MAX_ALT_COUNT;
import static com.hartwig.hmftools.sage.config.BaseQualityRecalibrationConfig.DEFAULT_BQR_SAMPLE_SIZE;
import static com.hartwig.hmftools.sage.config.FilterConfig.DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT;
import static com.hartwig.hmftools.sage.config.FilterConfig.DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY;
import static com.hartwig.hmftools.sage.config.FilterConfig.DEFAULT_HARD_MIN_TUMOR_QUAL;
import static com.hartwig.hmftools.sage.config.FilterConfig.DEFAULT_HIGH_CONFIDENCE_FILTER;
import static com.hartwig.hmftools.sage.config.FilterConfig.DEFAULT_HOTSPOT_FILTER;
import static com.hartwig.hmftools.sage.config.FilterConfig.DEFAULT_LOW_CONFIDENCE_FILTER;
import static com.hartwig.hmftools.sage.config.FilterConfig.DEFAULT_PANEL_FILTER;
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_BASE_QUAL_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_JITTER_MIN_REPEAT_COUNT;
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_JITTER_PENALTY;
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_MAP_QUAL_DISTANCE_FROM_REF;
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_MAP_QUAL_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY;
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_READ_EDGE_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_MAX_READ_DEPTH;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_MAX_READ_DEPTH_PANEL;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_MAX_REALIGNMENT_DEPTH;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_MIN_BASE_QUALITY;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_THREADS;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;

import org.jetbrains.annotations.NotNull;

public class SageConfigTest {

    @NotNull
    public static SageConfig testConfig() {
        return ImmutableSageConfig.builder()
                .panelOnly(false)
                .version("2.2")
                .outputFile("output file")
                .transcriptRegions(HmfGenePanelSupplier.allGeneList37())
                .reference(Lists.newArrayList("reference"))
                .referenceBam(Lists.newArrayList("referenceBam"))
                .tumor(Lists.newArrayList("tumorList"))
                .tumorBam(Lists.newArrayList("tumorBamList"))
                .refGenome("refGenome")
                .panelBed("panel")
                .highConfidenceBed("highConfidence")
                .hotspots("hotspots")
                .threads(DEFAULT_THREADS)
                .minMapQuality(DEFAULT_MIN_MAP_QUALITY)
                .minBaseQuality(DEFAULT_MIN_BASE_QUALITY)
                .maxRealignmentDepth(DEFAULT_MAX_REALIGNMENT_DEPTH)
                .maxReadDepth(DEFAULT_MAX_READ_DEPTH)
                .maxReadDepthPanel(DEFAULT_MAX_READ_DEPTH_PANEL)
                .qualityConfig(defaultQualityConfig())
                .filter(defaultFilterConfig())
                .baseQualityRecalibrationConfig(defaultQualityRecalibrationConfig())
                .build();

    }

    public static QualityConfig defaultQualityConfig() {
        return ImmutableQualityConfig.builder()
                .jitterPenalty(DEFAULT_JITTER_PENALTY)
                .jitterMinRepeatCount(DEFAULT_JITTER_MIN_REPEAT_COUNT)
                .baseQualityFixedPenalty(DEFAULT_BASE_QUAL_FIXED_PENALTY)
                .distanceFromReadEdgeFixedPenalty(DEFAULT_READ_EDGE_FIXED_PENALTY)
                .mapQualityFixedPenalty(DEFAULT_MAP_QUAL_FIXED_PENALTY)
                .mapQualityAdditionalDistanceFromRefPenalty(DEFAULT_MAP_QUAL_DISTANCE_FROM_REF)
                .mapQualityImproperPairPenalty(DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY)
                .build();
    }

    public static FilterConfig defaultFilterConfig() {
        return ImmutableFilterConfig.builder()
                .hardFilter(false)
                .disableSoftFilter(false)
                .hardMinTumorQual(DEFAULT_HARD_MIN_TUMOR_QUAL)
                .hardMinTumorRawAltSupport(DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT)
                .hardMinTumorRawBaseQuality(DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY)
                .softHotspotFilter(DEFAULT_HOTSPOT_FILTER)
                .softPanelFilter(DEFAULT_PANEL_FILTER)
                .softHighConfidenceFilter(DEFAULT_HIGH_CONFIDENCE_FILTER)
                .softLowConfidenceFilter(DEFAULT_LOW_CONFIDENCE_FILTER)
                .build();
    }

    public static BaseQualityRecalibrationConfig defaultQualityRecalibrationConfig() {
        return ImmutableBaseQualityRecalibrationConfig.builder()
                .enabled(false)
                .plot(false)
                .maxAltCount(DEFAULT_BQR_MAX_ALT_COUNT)
                .sampleSize(DEFAULT_BQR_SAMPLE_SIZE)
                .build();
    }

}
