package com.hartwig.hmftools.sage.common;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.internal.Sets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.SageConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.RegionResults;
import com.hartwig.hmftools.sage.pipeline.RegionTask;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;

public class RegionTaskTester
{
    public final RegionResults Results;
    public final SageConfig Config;
    public MockRefGenome RefGenome;

    public final List<VariantHotspot> Hotspots;
    public final List<BaseRegion> PanelRegions;
    public final List<TranscriptData> Transcripts;
    public final List<BaseRegion> HighConfidenceRegions;
    public final Map<String, QualityRecalibrationMap> QualityRecalibrationMap;
    public final PhaseSetCounter PhaseSetCounter;
    public final Coverage Coverage;
    public final SamSlicerFactory SamSlicerFactory;

    public final MockSamSlicer TumorSamSlicer;

    public static final String TEST_TUMOR_ID = "TUMOR_ID";
    public static final String TEST_REF_ID = "TEST_REF_ID";

    public RegionTaskTester()
    {
        Results = new RegionResults(null);
        Config = new SageConfig();
        RefGenome = new MockRefGenome();

        Hotspots = Lists.newArrayList();
        PanelRegions = Lists.newArrayList();
        Transcripts = Lists.newArrayList();
        HighConfidenceRegions = Lists.newArrayList();
        QualityRecalibrationMap = Maps.newHashMap();
        PhaseSetCounter = new PhaseSetCounter();
        Coverage = new Coverage(Sets.newHashSet(), Collections.EMPTY_LIST);

        SamSlicerFactory = new SamSlicerFactory();

        Config.TumorIds.add(TEST_TUMOR_ID);

        TumorSamSlicer = new MockSamSlicer();
        SamSlicerFactory.addSamSlicer(TEST_TUMOR_ID, TumorSamSlicer);
    }

    public RegionTask createRegionTask(final ChrBaseRegion region)
    {
        return new RegionTask(
                0, region, Results, Config, RefGenome, Hotspots, PanelRegions, Transcripts, HighConfidenceRegions,
                QualityRecalibrationMap, PhaseSetCounter, Coverage, SamSlicerFactory);
    }

}
