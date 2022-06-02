package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.common.TestUtils.createSamRecord;

import static junit.framework.TestCase.assertEquals;

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
import com.hartwig.hmftools.sage.common.MockSamSlicer;
import com.hartwig.hmftools.sage.common.SageVariant;
import com.hartwig.hmftools.sage.common.SamSlicerFactory;
import com.hartwig.hmftools.sage.coverage.Coverage;
import com.hartwig.hmftools.sage.phase.PhaseSetCounter;
import com.hartwig.hmftools.sage.pipeline.RegionResults;
import com.hartwig.hmftools.sage.pipeline.RegionTask;
import com.hartwig.hmftools.sage.quality.QualityRecalibrationMap;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class RegionTaskTest
{
    private final RegionResults mResults;
    private final SageConfig mConfig;
    private MockRefGenome mRefGenome;

    private final List<VariantHotspot> mHotspots;
    private final List<BaseRegion> mPanelRegions;
    private final List<TranscriptData> mTranscripts;
    private final List<BaseRegion> mHighConfidenceRegions;
    private final Map<String, QualityRecalibrationMap> mQualityRecalibrationMap;
    private final PhaseSetCounter mPhaseSetCounter;
    private final Coverage mCoverage;
    private final SamSlicerFactory mSamSlicerFactory;

    private final MockSamSlicer mTumorSamSlicer;

    private static final String TEST_TUMOR_ID = "TUMOR_ID";
    private static final String TEST_REF_ID = "TEST_REF_ID";

    private static final String TEST_REF_BASES =
            "GCAGGAGAATCCCTTGAACCTGGGAGGCAGAGGTTACAGTGAGCTGAGAT"
          + "CATGCCATTGCACTCTAGCCTGGGCAACAAGAGTGAAACTCCGCCTCAAA"
          + "ACAAACAAACAAACAAACAAACAAACAAACAAACAAAAACCTCCAAAACA";
           //0         10        20        30        40       49

    private final String mRefBases;

    public RegionTaskTest()
    {
        mResults = new RegionResults(null);
        mConfig = new SageConfig();
        mRefGenome = new MockRefGenome();

        mRefBases = TEST_REF_BASES; // generateRandomBases(200);
        mRefGenome.RefGenomeMap.put(CHR_1, "X" + mRefBases + generateRandomBases(1500)); // need to cover the ref sequence buffer

        mHotspots = Lists.newArrayList();
        mPanelRegions = Lists.newArrayList();
        mTranscripts = Lists.newArrayList();
        mHighConfidenceRegions = Lists.newArrayList();
        mQualityRecalibrationMap = Maps.newHashMap();
        mPhaseSetCounter = new PhaseSetCounter();
        mCoverage = new Coverage(Sets.newHashSet(), Collections.EMPTY_LIST);

        mSamSlicerFactory = new SamSlicerFactory();

        mConfig.TumorIds.add(TEST_TUMOR_ID);
        // mConfig.TumorBams.add("TumorBamFile");

        mTumorSamSlicer = new MockSamSlicer();
        mSamSlicerFactory.addSamSlicer(TEST_TUMOR_ID, mTumorSamSlicer);
    }

    @Test
    public void testBasicVariants()
    {
        ChrBaseRegion region = new ChrBaseRegion(CHR_1, 1, 200);
        RegionTask task = createRegionTask(region);

        String readBases = mRefBases.substring(0, 20) + "A" + mRefBases.substring(21, 51);

        SAMRecord read1 = createSamRecord("READ_01", CHR_1, 1, readBases, "50M");
        mTumorSamSlicer.ReadRecords.add(read1);

        readBases = mRefBases.substring(2, 20) + "A" + mRefBases.substring(21, 53);
        SAMRecord read2 = createSamRecord("READ_02", CHR_1, 3, readBases, "50M");
        mTumorSamSlicer.ReadRecords.add(read2);

        task.run();

        assertEquals(1, task.getVariants().size());
        SageVariant var = task.getVariants().get(0);
        assertEquals(21, var.position());
        assertEquals("T", var.ref());
        assertEquals("A", var.alt());
    }

    private RegionTask createRegionTask(final ChrBaseRegion region)
    {
        return new RegionTask(
                0, region, mResults, mConfig, mRefGenome, mHotspots, mPanelRegions, mTranscripts, mHighConfidenceRegions,
        mQualityRecalibrationMap, mPhaseSetCounter, mCoverage, mSamSlicerFactory);
    }

}
