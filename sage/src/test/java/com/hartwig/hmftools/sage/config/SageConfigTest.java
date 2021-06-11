package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.sage.config.BaseQualityRecalibrationConfig.DEFAULT_BQR_MAX_ALT_COUNT;
import static com.hartwig.hmftools.sage.config.BaseQualityRecalibrationConfig.DEFAULT_BQR_MIN_MAP_QUAL;
import static com.hartwig.hmftools.sage.config.BaseQualityRecalibrationConfig.DEFAULT_BQR_SAMPLE_SIZE;
import static com.hartwig.hmftools.sage.config.FilterConfig.DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT;
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
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_MAP_QUAL_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY;
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY;
import static com.hartwig.hmftools.sage.config.QualityConfig.DEFAULT_READ_EDGE_FIXED_PENALTY;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_MAX_READ_DEPTH;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_MAX_READ_DEPTH_PANEL;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_MAX_REALIGNMENT_DEPTH;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_MIN_MAP_QUALITY;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_MNV;
import static com.hartwig.hmftools.sage.config.SageConfig.DEFAULT_THREADS;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.ValidationStringency;

public class SageConfigTest
{

    @Test
    public void testBqrFile()
    {
        SageConfig config = testConfig();
        assertEquals("SAMPLE.sage.bqr.tsv", config.baseQualityRecalibrationFile("SAMPLE"));

        config = ImmutableSageConfig.builder().from(testConfig()).outputFile("./out.vcf").build();
        assertEquals("./SAMPLE.sage.bqr.tsv", config.baseQualityRecalibrationFile("SAMPLE"));
    }

    @NotNull
    public static SageConfig testConfig()
    {
        return ImmutableSageConfig.builder()
                .panelOnly(false)
                .version("2.2")
                .inputFile("in.vcf")
                .outputFile("out.vcf")
                .transcriptRegions(HmfGenePanelSupplier.allGeneList37())
                .reference(Lists.newArrayList("reference"))
                .referenceBam(Lists.newArrayList("referenceBam"))
                .tumor(Lists.newArrayList("tumorList"))
                .tumorBam(Lists.newArrayList("tumorBamList"))
                .refGenome("refGenome")
                .panelBed("panel")
                .highConfidenceBed("highConfidence")
                .hotspots("hotspots")
                .mnvEnabled(DEFAULT_MNV)
                .coverageBed("")
                .threads(DEFAULT_THREADS)
                .minMapQuality(DEFAULT_MIN_MAP_QUALITY)
                .maxRealignmentDepth(DEFAULT_MAX_REALIGNMENT_DEPTH)
                .maxReadDepth(DEFAULT_MAX_READ_DEPTH)
                .maxReadDepthPanel(DEFAULT_MAX_READ_DEPTH_PANEL)
                .qualityConfig(defaultQualityConfig())
                .regionSliceSize(500_000)
                .filter(defaultFilterConfig())
                .readContextFlankSize(SageConfig.DEFAULT_READ_CONTEXT_FLANK_SIZE)
                .baseQualityRecalibrationConfig(defaultQualityRecalibrationConfig())
                .validationStringency(ValidationStringency.DEFAULT_STRINGENCY)
                .build();
    }

    public static QualityConfig defaultQualityConfig()
    {
        return ImmutableQualityConfig.builder()
                .highlyPolymorphicGenes(Lists.newArrayList())
                .jitterPenalty(DEFAULT_JITTER_PENALTY)
                .jitterMinRepeatCount(DEFAULT_JITTER_MIN_REPEAT_COUNT)
                .baseQualityFixedPenalty(DEFAULT_BASE_QUAL_FIXED_PENALTY)
                .distanceFromReadEdgeFixedPenalty(DEFAULT_READ_EDGE_FIXED_PENALTY)
                .mapQualityFixedPenalty(DEFAULT_MAP_QUAL_FIXED_PENALTY)
                .mapQualityReadEventsPenalty(DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY)
                .mapQualityImproperPairPenalty(DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY)
                .build();
    }

    public static FilterConfig defaultFilterConfig()
    {
        return ImmutableFilterConfig.builder()
                .hardFilter(false)
                .softFilter(true)
                .mnvFilter(true)
                .hardMinTumorQual(DEFAULT_HARD_MIN_TUMOR_QUAL)
                .hardMinTumorRawAltSupport(DEFAULT_HARD_MIN_TUMOR_ALT_SUPPORT)
                .hardMinTumorRawBaseQuality(DEFAULT_HARD_MIN_TUMOR_BASE_QUALITY)
                .softHotspotFilter(DEFAULT_HOTSPOT_FILTER)
                .softPanelFilter(DEFAULT_PANEL_FILTER)
                .softHighConfidenceFilter(DEFAULT_HIGH_CONFIDENCE_FILTER)
                .softLowConfidenceFilter(DEFAULT_LOW_CONFIDENCE_FILTER)
                .filteredMaxNormalAltSupport(DEFAULT_FILTERED_MAX_NORMAL_ALT_SUPPORT)
                .build();
    }

    public static BaseQualityRecalibrationConfig defaultQualityRecalibrationConfig()
    {
        return ImmutableBaseQualityRecalibrationConfig.builder()
                .enabled(false)
                .plot(false)
                .maxAltCount(DEFAULT_BQR_MAX_ALT_COUNT)
                .sampleSize(DEFAULT_BQR_SAMPLE_SIZE)
                .minMapQuality(DEFAULT_BQR_MIN_MAP_QUAL)
                .build();
    }
}
