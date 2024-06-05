package com.hartwig.hmftools.sage.common;

import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.SageCallConfig;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.evidence.FragmentLengths;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.RegionResults;
import com.hartwig.hmftools.sage.pipeline.RegionTask;
import com.hartwig.hmftools.sage.bqr.BqrRecordMap;
import com.hartwig.hmftools.sage.quality.MsiJitterCalcs;

public class RegionTaskTester
{
    public final RegionResults Results;
    public final SageCallConfig Config;
    public MockRefGenome RefGenome;

    public final List<SimpleVariant> Hotspots;
    public final List<BaseRegion> PanelRegions;
    public final List<TranscriptData> Transcripts;
    public final List<BaseRegion> HighConfidenceRegions;
    public final Map<String, BqrRecordMap> QualityRecalibrationMap;
    public final MsiJitterCalcs JitterCalcs;
    public final PhaseSetCounter PhaseSetCounter;
    public final Coverage Coverage;
    public final SamSlicerFactory SamSlicerFactory;

    public final MockSamSlicer TumorSamSlicer;

    public static final String TEST_TUMOR_ID = "TUMOR_ID";
    public static final String TEST_REF_ID = "TEST_REF_ID";

    public RegionTaskTester()
    {
        Results = new RegionResults(null);
        Config = new SageCallConfig();
        RefGenome = new MockRefGenome();

        Hotspots = Lists.newArrayList();
        PanelRegions = Lists.newArrayList();
        Transcripts = Lists.newArrayList();
        HighConfidenceRegions = Lists.newArrayList();
        QualityRecalibrationMap = Maps.newHashMap();
        JitterCalcs = new MsiJitterCalcs();
        PhaseSetCounter = new PhaseSetCounter();
        Coverage = new Coverage(Lists.newArrayList(), Collections.EMPTY_LIST, Config.Common);

        SamSlicerFactory = new SamSlicerFactory();

        Config.TumorIds.add(TEST_TUMOR_ID);

        TumorSamSlicer = new MockSamSlicer();
        SamSlicerFactory.addSamSlicer(TEST_TUMOR_ID, TumorSamSlicer);
    }

    public RegionTask createRegionTask(final ChrBaseRegion region)
    {
        return new RegionTask(
                0, region, Results, Config, RefGenome, Hotspots, PanelRegions, Transcripts, HighConfidenceRegions,
                QualityRecalibrationMap, JitterCalcs, PhaseSetCounter, Coverage, SamSlicerFactory, new FragmentLengths(Config.Common));
    }
}
