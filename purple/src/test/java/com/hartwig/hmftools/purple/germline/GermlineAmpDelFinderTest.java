package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.common.purple.GermlineStatus.AMPLIFICATION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.SegmentSupport.BND;
import static com.hartwig.hmftools.common.purple.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS_FILTER;
import static com.hartwig.hmftools.purple.germline.GermlineAmpDelFinder.FILTER_CN_INCONSISTENCY;
import static com.hartwig.hmftools.purple.germline.GermlineAmpDelFinder.FILTER_COHORT_FREQ;
import static com.hartwig.hmftools.purple.germline.GermlineAmpDelFinder.FILTER_LOW_DEPTH_NO_SVS;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.mockito.ArgumentMatchers.any;
import static org.mockito.ArgumentMatchers.eq;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.driver.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.CommonVcfTags;
import com.hartwig.hmftools.purple.drivers.AmpDelRegionFrequency;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.junit.Before;
import org.junit.Test;
import org.mockito.Mockito;

public class GermlineAmpDelFinderTest
{
    private final GeneData geneData1 = geneData(CHR_1, "Id1_1", "GENE1_1", 1200, 1300, "K1_1");
    private final GeneData geneData2 = geneData(CHR_1, "Id1_2", "GENE1_2", 2200, 2300, "K1_2");
    private final GeneData geneData3 = geneData(CHR_2, "Id2_1", "GENE2_1", 1200, 1300, "K2_1");
    private final int transId1 = 1;
    private final ExonData exon1_1_1 = exonData(transId1, 1, 1210, 1218);
    private final ExonData exon1_1_2 = exonData(transId1, 2, 1220, 1228);
    private final TranscriptData td1_1 = transcriptData(geneData1.GeneId, transId1, true, exon1_1_1, exon1_1_2);
    private final int transId2 = 2;
    private final ExonData exon1_2_1 = exonData(transId2, 1, 2210, 2218);
    private final ExonData exon1_2_2 = exonData(transId2, 2, 2220, 2228);
    private final TranscriptData td1_2 = transcriptData(geneData2.GeneId, transId2, true, exon1_2_1, exon1_2_2);
    private final int transId3 = 3;
    private final ExonData exon2_1_1 = exonData(transId3, 1, 1210, 1218);
    private final ExonData exon2_1_2 = exonData(transId3, 2, 1220, 1228);
    private final TranscriptData td2_1 = transcriptData(geneData3.GeneId, transId3, true, exon2_1_1, exon2_1_2);
    private DriverGene mDriverGene1;
    private DriverGene mDriverGene2;
    private DriverGene driverGene3;

    private class GDS implements GermlineAmpDelFinder.GeneDataSupplier
    {
        List<TranscriptData> transcriptData = List.of(td1_1, td1_2, td2_1);
        Map<String, List<GeneData>> chrGeneMap = Map.of(CHR_1, List.of(geneData1, geneData2), CHR_2, List.of(geneData3));

        @Override
        public List<GeneData> getGeneData(final String chromosome)
        {
            return chrGeneMap.get(chromosome);
        }

        @Override
        public TranscriptData getTranscriptData(final String geneId)
        {
            return transcriptData.stream().filter(td -> td.GeneId.equals(geneId)).findFirst().orElse(null);
        }
    }

    private GDS mEnsemblDataCache;
    private GermlineAmpDelFrequencyCache mGermlineAmpDelFrequencyCache = Mockito.mock(GermlineAmpDelFrequencyCache.class);
    private GermlineAmpDelFinder mGermlineAmpDelFinder;

    @Before
    public void setup()
    {
        mDriverGene1 = driverGene(geneData1.GeneName, DriverGeneGermlineReporting.ANY);
        mDriverGene2 = driverGene(geneData2.GeneName, DriverGeneGermlineReporting.ANY);
        driverGene3 = driverGene(geneData3.GeneName, DriverGeneGermlineReporting.ANY);
        mEnsemblDataCache = new GDS();
        mGermlineAmpDelFrequencyCache = Mockito.mock(GermlineAmpDelFrequencyCache.class);
        Mockito.when(mGermlineAmpDelFrequencyCache.getRegionFrequency(
                any(String.class), any(Integer.class), any(Integer.class), any(Integer.class),
                any(AmpDelRegionFrequency.EventType.class))).thenReturn(3);
        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, mDriverGene1, geneData2.GeneName, mDriverGene2, geneData3.GeneName, driverGene3);
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);
    }

    @Test
    public void noStructuralVariants()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or = createObservedRegion(CHR_1, 1001, 1400, HOM_DELETION);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or), List.of());
        assertEquals(1, mGermlineAmpDelFinder.getEvents().size());
        assertEquals(mDriverGene1.gene(), mGermlineAmpDelFinder.getEvents().get(0).GeneName);
    }

    @Test
    public void copyNumbersDoNotIntersectGenes()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 5001, 6000, 2);
        ObservedRegion or = createObservedRegion(CHR_1, 1001, 1400, HOM_DELETION);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or), List.of());
        assertEquals(1, mGermlineAmpDelFinder.getEvents().size());
    }

    @Test
    public void observedRegionsDoNotIntersectGenes()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or = createObservedRegion(CHR_1, 6001, 7000, HOM_DELETION);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or), List.of());
        assertEquals(0, mGermlineAmpDelFinder.getEvents().size());
    }

    @Test
    public void multipleObservedRegions()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(CHR_1, 1001, 1400, HOM_DELETION);
        ObservedRegion or2 = createObservedRegion(CHR_1, 2101, 2400, HOM_DELETION);
        ObservedRegion or3 = createObservedRegion(CHR_2, 1001, 1400, HOM_DELETION);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1, or2, or3), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(3, deletions.size());
        assertEquals(mDriverGene1.gene(), deletions.get(0).GeneName);
        assertEquals(mDriverGene2.gene(), deletions.get(1).GeneName);
        assertEquals(driverGene3.gene(), deletions.get(2).GeneName);
    }

    @Test
    public void observedRegionCloseToDriverGene()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(CHR_1, 1001, 1100, HOM_DELETION);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        assertEquals(mDriverGene1.gene(), deletions.get(0).GeneName);
    }

    @Test
    public void singleDepthWindowAndNotBoundedByAStructuralVariant()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, NONE, NONE, 1, 2);
        ObservedRegion or1 = createObservedRegion(CHR_1, 1001, 1100, 1, HOM_DELETION);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> events = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, events.size());
        checkAllEventsPass(FILTER_LOW_DEPTH_NO_SVS, events, 0, ReportedStatus.NONE);
    }

    @Test
    public void multipleDepthWindowsAndNotBoundedByAStructuralVariant()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, NONE, NONE, 2, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1001, 1100, 2, HOM_DELETION, 0.2, DEFAULT_OBS_NORMAL_RATIO);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> events = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, events.size());
        checkAllEventsPass(PASS_FILTER, events, 0, ReportedStatus.REPORTED);
    }

    @Test
    public void singleDepthWindowAndBoundedByAStructuralVariantAtStart()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, BND, NONE, 1, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1001, 2000, 2, HOM_DELETION, 0.2, DEFAULT_OBS_NORMAL_RATIO);
        or1.setSupport(BND);
        ObservedRegion or2 = createObservedRegion(CHR_1, 2001, 3000, HOM_DELETION);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1, or2), List.of());
        List<GermlineAmpDel> events = mGermlineAmpDelFinder.getEvents();
        assertEquals(3, events.size());
        checkAllEventsPass(PASS_FILTER, events, 0, ReportedStatus.REPORTED);
    }

    @Test
    public void singleDepthWindowAndBoundedByAStructuralVariantAtEnd()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, NONE, BND, 1, 2);

        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1001, 2000, 2, HOM_DELETION, 0.2, DEFAULT_OBS_NORMAL_RATIO);

        ObservedRegion or2 = createObservedRegion(
                CHR_1, 2001, 3000, 2, HOM_DELETION, 0.2, DEFAULT_OBS_NORMAL_RATIO);

        or2.setSupport(BND);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1, or2), List.of());
        List<GermlineAmpDel> events = mGermlineAmpDelFinder.getEvents();
        assertEquals(3, events.size()); // both genes
        checkAllEventsPass(PASS_FILTER, events, 0, ReportedStatus.REPORTED);
    }

    @Test
    public void multipleDepthWindowsAndBoundedByAStructuralVariantAtStart()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, BND, NONE, 2, 0.2);

        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1001, 2000, 2, HOM_DELETION, 0.2, DEFAULT_OBS_NORMAL_RATIO);

        or1.setSupport(BND);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> events = mGermlineAmpDelFinder.getEvents();
        assertEquals(2, events.size()); // both genes
        checkAllEventsPass(PASS_FILTER, events, 0, ReportedStatus.REPORTED);
        checkAllEventsPass(PASS_FILTER, events, 1, ReportedStatus.REPORTED);
    }

    @Test
    public void changesThatAreNeitherAmpsNorDels()
    {
        PurpleCopyNumber pcn1 = createCopyNumber(CHR_1, 1001, 2000, 2);
        PurpleCopyNumber pcn2 = createCopyNumber(CHR_1, 2001, 3000, 2);
        ObservedRegion or1 = createObservedRegion(CHR_1, 1001, 1300, GermlineStatus.EXCLUDED);
        ObservedRegion or2 = createObservedRegion(CHR_1, 1301, 1600, GermlineStatus.DIPLOID);
        ObservedRegion or3 = createObservedRegion(CHR_1, 1601, 1900, GermlineStatus.CENTROMETIC);
        ObservedRegion or4 = createObservedRegion(CHR_1, 2001, 2300, GermlineStatus.UNKNOWN);
        ObservedRegion or5 = createObservedRegion(CHR_1, 2301, 2600, GermlineStatus.DIPLOID);
        ObservedRegion or6 = createObservedRegion(CHR_1, 2601, 2900, GermlineStatus.NOISE);
        mGermlineAmpDelFinder.findEvents(List.of(pcn1, pcn2), List.of(or1, or2, or3, or4, or5, or6), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(0, deletions.size());
    }

    @Test
    public void shortDeletion()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1201, 1500, HOM_DELETION, 1, 0.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        checkAllEventsPass(FILTER_LOW_DEPTH_NO_SVS, deletions, 0, ReportedStatus.NONE);
    }

    @Test
    public void shortAmp()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1201, 1500, AMPLIFICATION, 1, 1.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> events = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, events.size());
        checkAllEventsPass(FILTER_LOW_DEPTH_NO_SVS, events, 0, ReportedStatus.NONE);
    }

    private static void checkAllEventsPass(final String filterRegionLength, final List<GermlineAmpDel> events, final int index,
            final ReportedStatus none)
    {
        assertEquals(filterRegionLength, events.get(index).Filter);
        assertEquals(none, events.get(index).Reported);
    }

    @Test
    public void likelyDiploid()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, GermlineStatus.LIKELY_DIPLOID, 2, 0.8);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> events = mGermlineAmpDelFinder.getEvents();
        assertEquals(0, events.size());
    }

    @Test
    public void delThatIsFrequentInCohort()
    {
        mGermlineAmpDelFrequencyCache = Mockito.mock(GermlineAmpDelFrequencyCache.class);
        Mockito.when(mGermlineAmpDelFrequencyCache.getRegionFrequency(
                any(String.class), any(Integer.class), any(Integer.class), any(Integer.class),
                eq(AmpDelRegionFrequency.EventType.DEL))).thenReturn(4);
        Mockito.when(mGermlineAmpDelFrequencyCache.getRegionFrequency(
                any(String.class), any(Integer.class), any(Integer.class), any(Integer.class),
                eq(AmpDelRegionFrequency.EventType.AMP))).thenReturn(1);
        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, mDriverGene1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(CHR_1, List.of(geneData1));
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, HOM_DELETION, 1, 0.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        checkAllEventsPass(FILTER_COHORT_FREQ, deletions, 0, ReportedStatus.NONE);
    }

    @Test
    public void delThatIsInfrequentInCohort()
    {
        mGermlineAmpDelFrequencyCache = Mockito.mock(GermlineAmpDelFrequencyCache.class);
        Mockito.when(mGermlineAmpDelFrequencyCache.getRegionFrequency(
                any(String.class), any(Integer.class), any(Integer.class), any(Integer.class),
                eq(AmpDelRegionFrequency.EventType.DEL))).thenReturn(1);
        Mockito.when(mGermlineAmpDelFrequencyCache.getRegionFrequency(
                any(String.class), any(Integer.class), any(Integer.class), any(Integer.class),
                eq(AmpDelRegionFrequency.EventType.AMP))).thenReturn(10);
        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, mDriverGene1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(CHR_1, List.of(geneData1));
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, HOM_DELETION, 1, 0.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        checkAllEventsPass(PASS_FILTER, deletions, 0, ReportedStatus.REPORTED);
    }

    @Test
    public void ampThatIsFrequentInCohort()
    {
        mGermlineAmpDelFrequencyCache = Mockito.mock(GermlineAmpDelFrequencyCache.class);
        Mockito.when(mGermlineAmpDelFrequencyCache.getRegionFrequency(
                any(String.class), any(Integer.class), any(Integer.class), any(Integer.class),
                eq(AmpDelRegionFrequency.EventType.DEL))).thenReturn(1);
        Mockito.when(mGermlineAmpDelFrequencyCache.getRegionFrequency(
                any(String.class), any(Integer.class), any(Integer.class), any(Integer.class),
                eq(AmpDelRegionFrequency.EventType.AMP))).thenReturn(10);
        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, mDriverGene1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(CHR_1, List.of(geneData1));
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, AMPLIFICATION, 1, 1.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> events = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, events.size());
        checkAllEventsPass(FILTER_COHORT_FREQ, events, 0, ReportedStatus.NONE);
    }

    @Test
    public void ampThatIsInfrequentInCohort()
    {
        mGermlineAmpDelFrequencyCache = Mockito.mock(GermlineAmpDelFrequencyCache.class);
        Mockito.when(mGermlineAmpDelFrequencyCache.getRegionFrequency(
                any(String.class), any(Integer.class), any(Integer.class), any(Integer.class),
                eq(AmpDelRegionFrequency.EventType.DEL))).thenReturn(10);
        Mockito.when(mGermlineAmpDelFrequencyCache.getRegionFrequency(
                any(String.class), any(Integer.class), any(Integer.class), any(Integer.class),
                eq(AmpDelRegionFrequency.EventType.AMP))).thenReturn(1);
        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, mDriverGene1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(CHR_1, List.of(geneData1));
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, AMPLIFICATION, 1, 1.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        checkAllEventsPass(PASS_FILTER, deletions, 0, ReportedStatus.REPORTED);
    }

    @Test
    public void multipleFilters()
    {
        mGermlineAmpDelFrequencyCache = Mockito.mock(GermlineAmpDelFrequencyCache.class);
        Mockito.when(mGermlineAmpDelFrequencyCache.getRegionFrequency(
                any(String.class), any(Integer.class), any(Integer.class), any(Integer.class),
                eq(AmpDelRegionFrequency.EventType.DEL))).thenReturn(5);
        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, mDriverGene1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(CHR_1, List.of(geneData1));
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1201, 1500, 1, HOM_DELETION, 4, 0.25);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        String expectedFilter = String.format("%s;%s;%s", FILTER_LOW_DEPTH_NO_SVS, FILTER_CN_INCONSISTENCY, FILTER_COHORT_FREQ);
        checkAllEventsPass(expectedFilter, deletions, 0, ReportedStatus.NONE);
    }

    @Test
    public void copyNumberMajorAllele()
    {
        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, mDriverGene1, geneData3.GeneName, driverGene3);
        mEnsemblDataCache.transcriptData = List.of(td1_1, td2_1);
        mEnsemblDataCache.chrGeneMap = Map.of(CHR_1, List.of(geneData1), CHR_2, List.of(geneData3));
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn1 = createCopyNumber(CHR_1, 1001, 2000, 1);
        PurpleCopyNumber pcn2 = createCopyNumber(CHR_2, 1001, 2000, 0);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, HOM_DELETION, 1, 0.5);
        ObservedRegion or2 = createObservedRegion(
                CHR_2, 1000, 2001, 2, HOM_DELETION, 1, 0.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn1, pcn2), List.of(or1, or2), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(2, deletions.size());
        assertEquals(mDriverGene1.gene(), deletions.get(0).GeneName);
        checkAllEventsPass(CommonVcfTags.PASS_FILTER, deletions, 0, ReportedStatus.REPORTED);
        assertEquals(driverGene3.gene(), deletions.get(1).GeneName);
        checkAllEventsPass(FILTER_CN_INCONSISTENCY, deletions, 1, ReportedStatus.NONE);
    }

    @Test
    public void driverGeneReportVariantNotLost()
    {
        DriverGene driverGene = driverGene(geneData1.GeneName, DriverGeneGermlineReporting.VARIANT_NOT_LOST);
        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, driverGene);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(CHR_1, List.of(geneData1));
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, HOM_DELETION, 1, 0.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        checkAllEventsPass(CommonVcfTags.PASS_FILTER, deletions, 0, ReportedStatus.REPORTED);
    }

    @Test
    public void driverGeneReportNone()
    {
        DriverGene driverGene = driverGene(geneData1.GeneName, DriverGeneGermlineReporting.NONE);
        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, driverGene);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(CHR_1, List.of(geneData1));
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, HOM_DELETION, 1, 0.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        checkAllEventsPass(CommonVcfTags.PASS_FILTER, deletions, 0, ReportedStatus.NOT_REPORTED);
    }

    @Test
    public void driverGeneReportWildtypeLost()
    {
        DriverGene driverGene1 = driverGene(geneData1.GeneName, DriverGeneGermlineReporting.WILDTYPE_LOST);
        DriverGene driverGene2 = driverGene(geneData3.GeneName, DriverGeneGermlineReporting.WILDTYPE_LOST);

        Map<String, DriverGene> driverMap = Map.of(driverGene1.gene(), driverGene1, driverGene2.gene(), driverGene2);

        mEnsemblDataCache.transcriptData = List.of(td1_1, td2_1);
        mEnsemblDataCache.chrGeneMap = Map.of(CHR_1, List.of(geneData1), CHR_2, List.of(geneData3));
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn1 = createCopyNumber(CHR_1, 1001, 2000, 2);
        PurpleCopyNumber pcn2 = createCopyNumber(CHR_2, 1001, 2000, 1);

        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, HOM_DELETION, 1, 0.5);
        ObservedRegion or2 = createObservedRegion(
                CHR_2, 1000, 2001, 2, HET_DELETION, 0, 0.5);

        mGermlineAmpDelFinder.findEvents(List.of(pcn1, pcn2), List.of(or1, or2), List.of());

        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(2, deletions.size());
        assertEquals(driverGene1.gene(), deletions.get(0).GeneName);
        checkAllEventsPass(CommonVcfTags.PASS_FILTER, deletions, 0, ReportedStatus.REPORTED);
        assertEquals(driverGene3.gene(), deletions.get(1).GeneName);
        checkAllEventsPass(CommonVcfTags.PASS_FILTER, deletions, 1, ReportedStatus.NOT_REPORTED);
    }

    @Test
    public void reportedDeletion()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, HOM_DELETION, 1, 0.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(2, deletions.size());
        checkAllEventsPass(CommonVcfTags.PASS_FILTER, deletions, 0, ReportedStatus.REPORTED);
    }

    @Test
    public void reportGermlineAmplification()
    {
        DriverGene driverGene1 = ImmutableDriverGene.builder().from(mDriverGene1)
                .reportGermlineAmplification(false)
                .build();
        DriverGene driverGene2 = ImmutableDriverGene.builder().from(mDriverGene2)
                .reportGermlineAmplification(true)
                .build();

        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, driverGene1, geneData2.GeneName, driverGene2, geneData3.GeneName, driverGene3);
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, AMPLIFICATION, 1.4, 1.4);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> events = mGermlineAmpDelFinder.getEvents();
        assertEquals(2, events.size());
        assertEquals(driverGene1.gene(), events.get(0).GeneName);
        checkAllEventsPass(CommonVcfTags.PASS_FILTER, events, 0, ReportedStatus.NOT_REPORTED);
        assertEquals(driverGene2.gene(), events.get(1).GeneName);
        checkAllEventsPass(CommonVcfTags.PASS_FILTER, events, 1, ReportedStatus.REPORTED);
    }

    @Test
    public void singleDeletionStructuralVariant()
    {
        PurpleCopyNumber pcn0 = createCopyNumber(CHR_1, 1, 1000, 2);
        PurpleCopyNumber pcn1 = createCopyNumber(CHR_1, 1001, 2000, 2);
        PurpleCopyNumber pcn2 = createCopyNumber(CHR_1, 2001, 3000, 2);
        ObservedRegion or0 = createObservedRegion(
                CHR_1, 1, 1000, GermlineStatus.DIPLOID, 1, 0.5);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1001, 2000, HOM_DELETION, 1, 0.5);
        ObservedRegion or2 = createObservedRegion(
                CHR_1, 2001, 3000, GermlineStatus.DIPLOID, 1, 0.5);
        StructuralVariant sv =
                PurpleTestUtils.createStructuralVariant(CHR_1, 1201, CHR_1, 1300, StructuralVariantType.DEL, 0.5, 0.5).build();
        mGermlineAmpDelFinder.findEvents(List.of(pcn0, pcn1, pcn2), List.of(or0, or1, or2), List.of(sv));
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        assertEquals(sv.start().position(), deletions.get(0).RegionStart);
        assertEquals(Objects.requireNonNull(sv.end()).position(), deletions.get(0).RegionEnd);
    }

    @Test
    public void singleInversionStructuralVariant()
    {
        PurpleCopyNumber pcn0 = createCopyNumber(CHR_1, 1, 1000, 2);
        PurpleCopyNumber pcn1 = createCopyNumber(CHR_1, 1001, 2000, 2);
        PurpleCopyNumber pcn2 = createCopyNumber(CHR_1, 2001, 3000, 2);
        ObservedRegion or0 = createObservedRegion(
                CHR_1, 1, 1000, GermlineStatus.DIPLOID, 1, 0.5);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1001, 2000, HOM_DELETION, 1, 0.5);
        ObservedRegion or2 = createObservedRegion(
                CHR_1, 2001, 3000, GermlineStatus.DIPLOID, 1, 0.5);

        StructuralVariant sv =
                PurpleTestUtils.createStructuralVariant(CHR_1, 1201, CHR_1, 1300, StructuralVariantType.INV,
                        0.5, 0.5).build();
        mGermlineAmpDelFinder.findEvents(List.of(pcn0, pcn1, pcn2), List.of(or0, or1, or2), List.of(sv));
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(2, deletions.size());
        assertEquals(sv.start().position(), deletions.get(0).RegionStart);
        assertEquals(or1.end(), deletions.get(0).RegionEnd);
        assertEquals(geneData1.GeneName, deletions.get(0).GeneName);
        assertEquals(geneData1.KaryotypeBand, deletions.get(0).ChromosomeBand);
        assertEquals(sv.start().position(), deletions.get(1).RegionStart);
        assertEquals(or1.end(), deletions.get(1).RegionEnd);
        assertEquals(geneData2.GeneName, deletions.get(1).GeneName);
        assertEquals(geneData2.KaryotypeBand, deletions.get(1).ChromosomeBand);
    }

    @Test
    public void supStructuralVariant()
    {
        // Setup with only one gene to avoid multiple deletions
        Map<String, DriverGene> driverMap = Map.of(geneData1.GeneName, mDriverGene1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(CHR_1, List.of(geneData1));
        mGermlineAmpDelFinder = new GermlineAmpDelFinder(driverMap, mEnsemblDataCache, mGermlineAmpDelFrequencyCache);

        PurpleCopyNumber pcn0 = createCopyNumber(CHR_1, 1, 1000, 2);
        PurpleCopyNumber pcn1 = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or0 = createObservedRegion(
                CHR_1, 1, 1000, GermlineStatus.DIPLOID, 1, 0.5);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1001, 2000, AMPLIFICATION, 1, 0.5);

        // DUP
        StructuralVariant svDup =
                PurpleTestUtils.createStructuralVariant(CHR_1, 1201, CHR_1, 1300, StructuralVariantType.DUP, 0.5, 0.5).build();
        mGermlineAmpDelFinder.findEvents(List.of(pcn0, pcn1), List.of(or0, or1), List.of(svDup));
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        assertEquals(geneData1.GeneName, deletions.get(0).GeneName);
    }

    @Test
    public void multipleStructuralVariantsMatchingRegion()
    {
        PurpleCopyNumber pcn0 = createCopyNumber(CHR_1, 1, 1000, 2);
        PurpleCopyNumber pcn1 = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or0 = createObservedRegion(
                CHR_1, 1, 1000, GermlineStatus.DIPLOID, 1, 0.5);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1001, 2000, HOM_DELETION, 1, 0.5);

        // Create multiple SVs that match the region
        StructuralVariant svDel = PurpleTestUtils.createStructuralVariant(
                CHR_1, 1201, CHR_1, 1300, StructuralVariantType.DEL, 0.5, 0.5).build();
        StructuralVariant svInv = PurpleTestUtils.createStructuralVariant(
                CHR_1, 1250, CHR_1, 1350, StructuralVariantType.INV, 0.5, 0.5).build();
        StructuralVariant svDup = PurpleTestUtils.createStructuralVariant(
                CHR_1, 1150, CHR_1, 1400, StructuralVariantType.DUP, 0.5, 0.5).build();

        // DEL type should be preferred over other types
        mGermlineAmpDelFinder.findEvents(List.of(pcn0, pcn1), List.of(or0, or1), List.of(svInv, svDel, svDup));
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(1, deletions.size());
        assertEquals(svDel.start().position(), deletions.get(0).RegionStart);
        assertEquals(svDel.end().position(), deletions.get(0).RegionEnd);
    }

    @Test
    public void useSVsForBoundaries()
    {
        PurpleCopyNumber pcn0 = createCopyNumber(CHR_1, 1, 1000, 2);
        PurpleCopyNumber pcn1 = createCopyNumber(CHR_1, 1001, 2000, 2);
        PurpleCopyNumber pcn2 = createCopyNumber(CHR_1, 2001, 3000, 2);

        ObservedRegion or0 = createObservedRegion(
                CHR_1, 1, 1000, GermlineStatus.DIPLOID, 1, 0.5);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1001, 2000, HOM_DELETION, 1, 0.5);
        ObservedRegion or2 = createObservedRegion(
                CHR_1, 2001, 3000, HOM_DELETION, 1, 0.5);

        StructuralVariant sv1 =
                PurpleTestUtils.createStructuralVariant(CHR_1, exon1_1_1.Start + 2, CHR_1, exon1_1_1.End, StructuralVariantType.DEL, 0.5, 0.5)
                        .build();
        StructuralVariant sv2 =
                PurpleTestUtils.createStructuralVariant(CHR_1,
                        exon1_2_1.Start - 20, CHR_1, exon1_2_1.Start + 2, StructuralVariantType.DEL, 0.5, 0.5).build();

        mGermlineAmpDelFinder.findEvents(List.of(pcn0, pcn1, pcn2), List.of(or0, or1, or2), List.of(sv1, sv2));
        List<GermlineAmpDel> events = mGermlineAmpDelFinder.getEvents();
        assertEquals(3, events.size());
        assertEquals(sv1.start().position(), events.get(0).RegionStart);
        assertEquals(sv1.end().position(), events.get(0).RegionEnd);
    }

    @Test
    public void hetDeletionGermlineStatus()
    {
        PurpleCopyNumber pcn = createCopyNumber(CHR_1, 1001, 2000, 2);
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, HET_DELETION, 1, 0.5);
        mGermlineAmpDelFinder.findEvents(List.of(pcn), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(2, deletions.size());
        checkAllEventsPass(CommonVcfTags.PASS_FILTER, deletions, 0, ReportedStatus.REPORTED);
        assertEquals(HET_DELETION, deletions.get(0).NormalStatus);
    }

    @Test
    public void emptyCopyNumbersList()
    {
        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, 2, HOM_DELETION, 1, 0.5);
        mGermlineAmpDelFinder.findEvents(List.of(), List.of(or1), List.of());
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(2, deletions.size());
        checkAllEventsPass(CommonVcfTags.PASS_FILTER, deletions, 0, ReportedStatus.REPORTED);
    }

    @Test
    public void multipleChromosomesWithStructuralVariants()
    {
        PurpleCopyNumber pcn1 = createCopyNumber(CHR_1, 1001, 2000, 2);
        PurpleCopyNumber pcn2 = createCopyNumber(CHR_2, 1001, 2000, 2);

        ObservedRegion or1 = createObservedRegion(
                CHR_1, 1000, 2001, HOM_DELETION, 1, 0.5);
        ObservedRegion or2 = createObservedRegion(
                CHR_2, 1000, 2001, HOM_DELETION, 1, 0.5);

        StructuralVariant sv1 = PurpleTestUtils.createStructuralVariant(
                        CHR_1, 1201, CHR_1, 1300, StructuralVariantType.DEL, 0.5, 0.5).build();
        StructuralVariant sv2 = PurpleTestUtils.createStructuralVariant(
                CHR_2, 1201, CHR_2, 1300, StructuralVariantType.DEL, 0.5, 0.5).build();

        mGermlineAmpDelFinder.findEvents(List.of(pcn1, pcn2), List.of(or1, or2), List.of(sv1, sv2));
        List<GermlineAmpDel> deletions = mGermlineAmpDelFinder.getEvents();
        assertEquals(2, deletions.size());

        // Check that we have one deletion from each chromosome
        boolean foundChr1 = false;
        boolean foundChr2 = false;
        for(GermlineAmpDel deletion : deletions)
        {
            if(deletion.Chromosome.equals(CHR_1))
            {
                foundChr1 = true;
                assertEquals(sv1.start().position(), deletion.RegionStart);
                assertEquals(sv1.end().position(), deletion.RegionEnd);
            }
            else if(deletion.Chromosome.equals(CHR_2))
            {
                foundChr2 = true;
                assertEquals(sv2.start().position(), deletion.RegionStart);
                assertEquals(sv2.end().position(), deletion.RegionEnd);
            }
        }

        assertTrue(foundChr1);
        assertTrue(foundChr2);
    }

    private PurpleCopyNumber createCopyNumber(final String chromosome, int start, int end,
            SegmentSupport startSupport,
            SegmentSupport endSupport,
            int depthWindowCount,
            double majorAlleleCopyNumber)
    {
        ImmutablePurpleCopyNumber result = ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(start)
                .end(end)
                .averageTumorCopyNumber(2 * majorAlleleCopyNumber)
                .segmentStartSupport(startSupport)
                .segmentEndSupport(endSupport)
                .method(CopyNumberMethod.UNKNOWN)
                .bafCount(0)
                .depthWindowCount(depthWindowCount)
                .gcContent(0)
                .minStart(start)
                .maxStart(start)
                .averageObservedBAF(0.5)
                .averageActualBAF(0.5)
                .build();
        Preconditions.checkState(result.majorAlleleCopyNumber() == majorAlleleCopyNumber);
        return result;
    }

    private PurpleCopyNumber createCopyNumber(final String chromosome, int start, int end, double majorAlleleCopyNumber)
    {
        ImmutablePurpleCopyNumber result =
                PurpleTestUtils.createCopyNumber(chromosome, start, end, 2 * majorAlleleCopyNumber).build();
        Preconditions.checkState(result.majorAlleleCopyNumber() == majorAlleleCopyNumber);
        return result;
    }

    private static final double DEFAULT_OBS_NORMAL_RATIO = 0.9;
    private static final double DEFAULT_COPY_NUMBER = 2.0;

    private ObservedRegion createObservedRegion(final String chromosome, int start, int end, int depthWindowCount,
            final GermlineStatus status)
    {
        return createObservedRegion(
                chromosome, start, end, depthWindowCount, status, DEFAULT_COPY_NUMBER, DEFAULT_OBS_NORMAL_RATIO);
    }

    private ObservedRegion createObservedRegion(final String chromosome, int start, int end, final GermlineStatus status)
    {
        return createObservedRegion(chromosome, start, end, status, DEFAULT_COPY_NUMBER, DEFAULT_OBS_NORMAL_RATIO);
    }

    private ObservedRegion createObservedRegion(
            final String chromosome, int start, int end, final GermlineStatus status, double refNormalisedCopyNumber,
            double observedNormalRatio)
    {
        return createObservedRegion(chromosome, start, end, 1, status, refNormalisedCopyNumber, observedNormalRatio);
    }

    private ObservedRegion createObservedRegion(
            final String chromosome, int start, int end, int depthWindowCount, final GermlineStatus status,
            double refNormalisedCopyNumber, double observedNormalRatio)
    {
        return new ObservedRegion(
                chromosome,
                start,
                end,
                false,
                null, // SegmentSupport
                0,   // BAF count
                0.0, // observed BAF
                depthWindowCount,   // depth window count
                0.0, // observed tumor ratio
                observedNormalRatio,
                0.0, // unnormalisedObservedNormalRatio
                status,
                false, // svCluster
                0.50, // gc content
                start,
                start,
                0.0,
                0.0,
                0.0, // deviationPenalty
                0.0, // eventPenalty
                refNormalisedCopyNumber,
                0.0,
                0.0,
                0.0,
                0.0);
    }

    private GeneData geneData(String chromosome, String id, String name, int start, int end, String karyotypeBand)
    {
        GeneData result = new GeneData(id, name, chromosome, (byte) 0, start, end, karyotypeBand);
        Preconditions.checkState(Objects.equals(result.GeneId, id));
        Preconditions.checkState(Objects.equals(result.GeneName, name));
        Preconditions.checkState(Objects.equals(result.Chromosome, chromosome));
        Preconditions.checkState(result.GeneStart == start);
        Preconditions.checkState(result.GeneEnd == end);
        Preconditions.checkState(Objects.equals(result.KaryotypeBand, karyotypeBand));
        return result;
    }

    private ExonData exonData(int transcriptId, int exonNumber, int start, int end)
    {
        ExonData result = new ExonData(transcriptId, start, end, exonNumber, 0, 0);
        Preconditions.checkState(result.TransId == transcriptId);
        Preconditions.checkState(result.Start == start);
        Preconditions.checkState(result.End == end);
        Preconditions.checkState(result.Rank == exonNumber);
        return result;
    }

    private TranscriptData transcriptData(String geneId, int transcriptId, boolean isCanonical, ExonData... exons)
    {
        TranscriptData result = new TranscriptData(
                transcriptId, "TD" + transcriptId, geneId, isCanonical,
                (byte) 0, 0, 0, 0, 0, null, null);
        Preconditions.checkState(Objects.equals(result.GeneId, geneId));
        Preconditions.checkState(result.TransId == transcriptId);
        Preconditions.checkState(result.IsCanonical == isCanonical);
        List<ExonData> exonDataList = List.of(exons);
        exonDataList.forEach(exonData -> Preconditions.checkState(Objects.equals(exonData.TransId, transcriptId)));
        result.setExons(exonDataList);
        return result;
    }

    private DriverGene driverGene(String geneName, DriverGeneGermlineReporting reporting)
    {
        ImmutableDriverGene result = ImmutableDriverGene.builder()
                .gene(geneName)
                .reportGermlineDeletion(reporting)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportDeletion(true)
                .reportGermlineVariant(DriverGeneGermlineReporting.ANY)
                .reportGermlineHotspot(DriverGeneGermlineReporting.ANY)
                .reportSplice(true)
                .reportHetDeletion(true)
                .hetDeletionThreshold(0.23)
                .reportDisruption(true)
                .reportAmplification(true)
                .reportLoh(true)
                .reportGermlineAmplification(true)
                .amplificationRatio(0.23)
                .reportSomaticHotspot(true)
                .likelihoodType(DriverCategory.TSG)
                .reportGermlineDisruption(DriverGeneGermlineReporting.ANY)
                .reportPGX(true)
                .reportHighExpression(false)
                .reportLowExpression(false)
                .reportNovelSpliceJunction(false)
                .build();
        Preconditions.checkState(Objects.equals(result.gene(), geneName));
        Preconditions.checkState(Objects.equals(result.reportGermlineDeletion(), reporting));
        return result;
    }
}
