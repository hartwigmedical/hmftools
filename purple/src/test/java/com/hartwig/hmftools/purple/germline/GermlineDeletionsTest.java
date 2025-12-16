package com.hartwig.hmftools.purple.germline;

import static com.hartwig.hmftools.purple.germline.GermlineDeletions.FILTER_CN_INCONSISTENCY;
import static com.hartwig.hmftools.purple.germline.GermlineDeletions.FILTER_COHORT_FREQ;
import static com.hartwig.hmftools.purple.germline.GermlineDeletions.FILTER_REGION_LENGTH;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.mockito.ArgumentMatchers.any;

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
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.CommonVcfTags;
import com.hartwig.hmftools.purple.region.ObservedRegion;

import org.junit.Before;
import org.junit.Test;
import org.mockito.Mockito;

public class GermlineDeletionsTest
{
    private final String chr1 = "chr1";
    private final String chr2 = "chr2";
    private final GeneData gd1_1 = geneData(chr1, "Id1_1", "GENE1_1", 1200, 1300, "K1_1");
    private final GeneData gd1_2 = geneData(chr1, "Id1_2", "GENE1_2", 2200, 2300, "K1_2");
    private final GeneData gd2_1 = geneData(chr2, "Id2_1", "GENE2_1", 1200, 1300, "K2_1");
    private final int transcriptId1_1 = 1;
    private final ExonData exon1_1_1 = exonData(transcriptId1_1, 1, 1210, 1218);
    private final ExonData exon1_1_2 = exonData(transcriptId1_1, 2, 1220, 1228);
    private final TranscriptData td1_1 = transcriptData(gd1_1.GeneId, transcriptId1_1, true, exon1_1_1, exon1_1_2);
    private final int transcriptId1_2 = 2;
    private final ExonData exon1_2_1 = exonData(transcriptId1_2, 1, 2210, 2218);
    private final ExonData exon1_2_2 = exonData(transcriptId1_2, 2, 2220, 2228);
    private final TranscriptData td1_2 = transcriptData(gd1_2.GeneId, transcriptId1_2, true, exon1_2_1, exon1_2_2);
    private final int transcriptId2_1 = 3;
    private final ExonData exon2_1_1 = exonData(transcriptId2_1, 1, 1210, 1218);
    private final ExonData exon2_1_2 = exonData(transcriptId2_1, 2, 1220, 1228);
    private final TranscriptData td2_1 = transcriptData(gd2_1.GeneId, transcriptId2_1, true, exon2_1_1, exon2_1_2);
    private final DriverGene driver1_1 = driverGene(gd1_1.GeneName, DriverGeneGermlineReporting.ANY);
    private final DriverGene driver1_2 = driverGene(gd1_2.GeneName, DriverGeneGermlineReporting.ANY);
    private final DriverGene driver2_1 = driverGene(gd2_1.GeneName, DriverGeneGermlineReporting.ANY);

    private class GDS implements GermlineDeletions.GeneDataSupplier
    {
        List<TranscriptData> transcriptData = List.of(td1_1, td1_2, td2_1);
        Map<String, List<GeneData>> chrGeneMap = Map.of(chr1, List.of(gd1_1, gd1_2), chr2, List.of(gd2_1));
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
    private GermlineDeletionFrequency mGermlineDeletionFrequency = Mockito.mock(GermlineDeletionFrequency.class);
    private GermlineDeletions mGermlineDeletions;

    @Before
    public void setup()
    {
        mEnsemblDataCache = new GDS();
        mGermlineDeletionFrequency = Mockito.mock(GermlineDeletionFrequency.class);
        Mockito.when(mGermlineDeletionFrequency.getRegionFrequency(any(String.class), any(Integer.class), any(Integer.class), any(Integer.class))).thenReturn(3);
        Map<String,DriverGene> driverMap = Map.of(gd1_1.GeneName, driver1_1, gd1_2.GeneName, driver1_2, gd2_1.GeneName, driver2_1);
        mGermlineDeletions = new GermlineDeletions(driverMap, mEnsemblDataCache, mGermlineDeletionFrequency);
    }

    @Test
    public void noStructuralVariants()
    {
        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or = or(chr1, 1001, 1400, GermlineStatus.HOM_DELETION);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or), List.of());
        assertEquals(1, mGermlineDeletions.getDeletions().size());
        assertEquals(driver1_1.gene(), mGermlineDeletions.getDeletions().get(0).GeneName);
    }

    @Test
    public void copyNumbersDoNotIntersectGenes()
    {
        final PurpleCopyNumber pcn = pcn(chr1, 5001, 6000, 2);
        final ObservedRegion or = or(chr1, 1001, 1400, GermlineStatus.HOM_DELETION);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or), List.of());
        assertEquals(1, mGermlineDeletions.getDeletions().size());
    }

    @Test
    public void observedRegionsDoNotIntersectGenes()
    {
        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or = or(chr1, 6001, 7000, GermlineStatus.HOM_DELETION);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or), List.of());
        assertEquals(0, mGermlineDeletions.getDeletions().size());
    }

    @Test
    public void multipleObservedRegions()
    {
        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1001, 1400, GermlineStatus.HOM_DELETION);
        final ObservedRegion or2 = or(chr1, 2101, 2400, GermlineStatus.HOM_DELETION);
        final ObservedRegion or3 = or(chr2, 1001, 1400, GermlineStatus.HOM_DELETION);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or1, or2, or3), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(3, deletions.size());
        assertEquals(driver1_1.gene(), deletions.get(0).GeneName);
        assertEquals(driver1_2.gene(), deletions.get(1).GeneName);
        assertEquals(driver2_1.gene(), deletions.get(2).GeneName);
    }

    @Test
    public void observedRegionCloseToDriverGene()
    {
        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1001, 1100, GermlineStatus.HOM_DELETION);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or1), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        assertEquals(driver1_1.gene(), deletions.get(0).GeneName);
    }

    @Test
    public void nonDeletions()
    {
        final PurpleCopyNumber pcn1 = pcn(chr1, 1001, 2000, 2);
        final PurpleCopyNumber pcn2 = pcn(chr1, 2001, 3000, 2);
        final ObservedRegion or1 = or(chr1, 1001, 1300, GermlineStatus.EXCLUDED);
        final ObservedRegion or2 = or(chr1, 1301, 1600, GermlineStatus.DIPLOID);
        final ObservedRegion or3 = or(chr1, 1601, 1900, GermlineStatus.CENTROMETIC);
        final ObservedRegion or4 = or(chr1, 2001, 2300, GermlineStatus.UNKNOWN);
        final ObservedRegion or5 = or(chr1, 2301, 2600, GermlineStatus.AMPLIFICATION);
        final ObservedRegion or6 = or(chr1, 2601, 2900, GermlineStatus.NOISE);
        mGermlineDeletions.findDeletions(List.of(pcn1, pcn2), List.of(or1, or2, or3,or4, or5, or6), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(0, deletions.size());
    }

    @Test
    public void shortDeletion()
    {
        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1201, 1500, GermlineStatus.HOM_DELETION, 1201, 1201, 1, 0.5);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or1), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        assertEquals(FILTER_REGION_LENGTH, deletions.get(0).Filter);
        assertEquals(ReportedStatus.NONE, deletions.get(0).Reported);
    }

    @Test
    public void inconsistentCopyNumber()
    {
        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 2, 0.9);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or1), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(2, deletions.size());
        assertEquals(FILTER_CN_INCONSISTENCY, deletions.get(0).Filter);
        assertEquals(ReportedStatus.NONE, deletions.get(0).Reported);
    }

    @Test
    public void cohortFrequency()
    {
        mGermlineDeletionFrequency = Mockito.mock(GermlineDeletionFrequency.class);
        Mockito.when(mGermlineDeletionFrequency.getRegionFrequency(any(String.class), any(Integer.class), any(Integer.class), any(Integer.class))).thenReturn(4);
        Map<String,DriverGene> driverMap = Map.of(gd1_1.GeneName, driver1_1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(chr1, List.of(gd1_1));
        mGermlineDeletions = new GermlineDeletions(driverMap, mEnsemblDataCache, mGermlineDeletionFrequency);

        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 1, 0.5);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or1), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        assertEquals(FILTER_COHORT_FREQ, deletions.get(0).Filter);
        assertEquals(ReportedStatus.NONE, deletions.get(0).Reported);
    }

    @Test
    public void multipleFilters()
    {
        mGermlineDeletionFrequency = Mockito.mock(GermlineDeletionFrequency.class);
        Mockito.when(mGermlineDeletionFrequency.getRegionFrequency(any(String.class), any(Integer.class), any(Integer.class), any(Integer.class))).thenReturn(5);
        Map<String,DriverGene> driverMap = Map.of(gd1_1.GeneName, driver1_1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(chr1, List.of(gd1_1));
        mGermlineDeletions = new GermlineDeletions(driverMap, mEnsemblDataCache, mGermlineDeletionFrequency);

        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1201, 1500, GermlineStatus.HOM_DELETION, 1201, 1201, 4, 0.5);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or1), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        String expectedFilter = String.format("%s;%s;%s", FILTER_REGION_LENGTH, FILTER_CN_INCONSISTENCY, FILTER_COHORT_FREQ);
        assertEquals(expectedFilter, deletions.get(0).Filter);
        assertEquals(ReportedStatus.NONE, deletions.get(0).Reported);
    }

    @Test
    public void copyNumberMajorAllele()
    {
        Map<String,DriverGene> driverMap = Map.of(gd1_1.GeneName, driver1_1, gd2_1.GeneName, driver2_1);
        mEnsemblDataCache.transcriptData = List.of(td1_1, td2_1);
        mEnsemblDataCache.chrGeneMap = Map.of(chr1, List.of(gd1_1), chr2, List.of(gd2_1));
        mGermlineDeletions = new GermlineDeletions(driverMap, mEnsemblDataCache, mGermlineDeletionFrequency);

        final PurpleCopyNumber pcn1 = pcn(chr1, 1001, 2000, 1);
        final PurpleCopyNumber pcn2 = pcn(chr2, 1001, 2000, 0);
        final ObservedRegion or1 = or(chr1, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 1, 0.5);
        final ObservedRegion or2 = or(chr2, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 1, 0.5);
        mGermlineDeletions.findDeletions(List.of(pcn1, pcn2), List.of(or1, or2), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(2, deletions.size());
        assertEquals(driver1_1.gene(), deletions.get(0).GeneName);
        assertEquals(CommonVcfTags.PASS, deletions.get(0).Filter);
        assertEquals(ReportedStatus.REPORTED, deletions.get(0).Reported);
        assertEquals(driver2_1.gene(), deletions.get(1).GeneName);
        assertEquals(FILTER_CN_INCONSISTENCY, deletions.get(1).Filter);
        assertEquals(ReportedStatus.NONE, deletions.get(1).Reported);
    }

    @Test
    public void driverGeneReportVariantNotLost()
    {
        DriverGene driverGene = driverGene(gd1_1.GeneName, DriverGeneGermlineReporting.VARIANT_NOT_LOST);
        Map<String, DriverGene> driverMap = Map.of(gd1_1.GeneName, driverGene);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(chr1, List.of(gd1_1));
        mGermlineDeletions = new GermlineDeletions(driverMap, mEnsemblDataCache, mGermlineDeletionFrequency);

        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 1, 0.5);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or1), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        assertEquals(CommonVcfTags.PASS, deletions.get(0).Filter);
        assertEquals(ReportedStatus.REPORTED, deletions.get(0).Reported);
    }

    @Test
    public void driverGeneReportNone()
    {
        DriverGene driverGene = driverGene(gd1_1.GeneName, DriverGeneGermlineReporting.NONE);
        Map<String, DriverGene> driverMap = Map.of(gd1_1.GeneName, driverGene);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(chr1, List.of(gd1_1));
        mGermlineDeletions = new GermlineDeletions(driverMap, mEnsemblDataCache, mGermlineDeletionFrequency);

        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 1, 0.5);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or1), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        assertEquals(CommonVcfTags.PASS, deletions.get(0).Filter);
        assertEquals(ReportedStatus.NONE, deletions.get(0).Reported);
    }

    @Test
    public void driverGeneReportWildtypeLost()
    {
        DriverGene driverGene1 = driverGene(gd1_1.GeneName, DriverGeneGermlineReporting.WILDTYPE_LOST);
        DriverGene driverGene2 = driverGene(gd2_1.GeneName, DriverGeneGermlineReporting.WILDTYPE_LOST);
        Map<String, DriverGene> driverMap = Map.of(gd1_1.GeneName, driverGene1, gd2_1.GeneName, driverGene2);
        mEnsemblDataCache.transcriptData = List.of(td1_1, td2_1);
        mEnsemblDataCache.chrGeneMap = Map.of(chr1, List.of(gd1_1), chr2, List.of(gd2_1));
        mGermlineDeletions = new GermlineDeletions(driverMap, mEnsemblDataCache, mGermlineDeletionFrequency);

        final PurpleCopyNumber pcn1 = pcn(chr1, 1001, 2000, 2);
        final PurpleCopyNumber pcn2 = pcn(chr2, 1001, 2000, 1);
        final ObservedRegion or1 = or(chr1, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 1, 0.5);
        final ObservedRegion or2 = or(chr2, 1000, 2001, GermlineStatus.HET_DELETION, 1000, 1000, 0, 0.5);
        mGermlineDeletions.findDeletions(List.of(pcn1, pcn2), List.of(or1, or2), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(2, deletions.size());
        assertEquals(driver1_1.gene(), deletions.get(0).GeneName);
        assertEquals(CommonVcfTags.PASS, deletions.get(0).Filter);
        assertEquals(ReportedStatus.NONE, deletions.get(0).Reported);
        assertEquals(driver2_1.gene(), deletions.get(1).GeneName);
        assertEquals(CommonVcfTags.PASS, deletions.get(1).Filter);
        assertEquals(ReportedStatus.REPORTED, deletions.get(1).Reported);
    }

    @Test
    public void reportedDeletion()
    {
        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 1, 0.5);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or1), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(2, deletions.size());
        assertEquals(CommonVcfTags.PASS, deletions.get(0).Filter);
        assertEquals(ReportedStatus.REPORTED, deletions.get(0).Reported);
    }

    @Test
    public void singleDeletionStructuralVariant()
    {
        final PurpleCopyNumber pcn0 = pcn(chr1, 1, 1000, 2);
        final PurpleCopyNumber pcn1 = pcn(chr1, 1001, 2000, 2);
        final PurpleCopyNumber pcn2 = pcn(chr1, 2001, 3000, 2);
        final ObservedRegion or0 = or(chr1, 1, 1000, GermlineStatus.DIPLOID, 1, 1, 1, 0.5);
        final ObservedRegion or1 = or(chr1, 1001, 2000, GermlineStatus.HOM_DELETION, 1001, 1001, 1, 0.5);
        final ObservedRegion or2 = or(chr1, 2001, 3000, GermlineStatus.DIPLOID, 2001, 2001, 1, 0.5);
        final StructuralVariant sv = PurpleTestUtils.createStructuralVariant(chr1, 1201, chr1, 1300, StructuralVariantType.DEL, 0.5, 0.5).build();
        mGermlineDeletions.findDeletions(List.of(pcn0, pcn1, pcn2), List.of(or0,or1,or2), List.of(sv));
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        assertEquals(sv.start().position(), deletions.get(0).RegionStart);
        assertEquals(Objects.requireNonNull(sv.end()).position(), deletions.get(0).RegionEnd);
    }

    @Test
    public void singleInversionStructuralVariant()
    {
        final PurpleCopyNumber pcn0 = pcn(chr1, 1, 1000, 2);
        final PurpleCopyNumber pcn1 = pcn(chr1, 1001, 2000, 2);
        final PurpleCopyNumber pcn2 = pcn(chr1, 2001, 3000, 2);
        final ObservedRegion or0 = or(chr1, 1, 1000, GermlineStatus.DIPLOID, 1, 1, 1, 0.5);
        final ObservedRegion or1 = or(chr1, 1001, 2000, GermlineStatus.HOM_DELETION, 1001, 1001, 1, 0.5);
        final ObservedRegion or2 = or(chr1, 2001, 3000, GermlineStatus.DIPLOID, 2001, 2001, 1, 0.5);
        final StructuralVariant sv = PurpleTestUtils.createStructuralVariant(chr1, 1201, chr1, 1300, StructuralVariantType.INV, 0.5, 0.5).build();
        mGermlineDeletions.findDeletions(List.of(pcn0, pcn1, pcn2), List.of(or0,or1,or2), List.of(sv));
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(2, deletions.size());
        assertEquals(sv.start().position(), deletions.get(0).RegionStart);
        assertEquals(or1.end(), deletions.get(0).RegionEnd);
        assertEquals(gd1_1.GeneName, deletions.get(0).GeneName);
        assertEquals(gd1_1.KaryotypeBand, deletions.get(0).ChromosomeBand);
        assertEquals(sv.start().position(), deletions.get(1).RegionStart);
        assertEquals(or1.end(), deletions.get(1).RegionEnd);
        assertEquals(gd1_2.GeneName, deletions.get(1).GeneName);
        assertEquals(gd1_2.KaryotypeBand, deletions.get(1).ChromosomeBand);
    }

    @Test
    public void otherStructuralVariantTypes()
    {
        // Setup with only one gene to avoid multiple deletions
        Map<String,DriverGene> driverMap = Map.of(gd1_1.GeneName, driver1_1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(chr1, List.of(gd1_1));
        mGermlineDeletions = new GermlineDeletions(driverMap, mEnsemblDataCache, mGermlineDeletionFrequency);

        final PurpleCopyNumber pcn0 = pcn(chr1, 1, 1000, 2);
        final PurpleCopyNumber pcn1 = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or0 = or(chr1, 1, 1000, GermlineStatus.DIPLOID, 1, 1, 1, 0.5);
        final ObservedRegion or1 = or(chr1, 1001, 2000, GermlineStatus.HOM_DELETION, 1001, 1001, 1, 0.5);

        // DUP
        final StructuralVariant svDup = PurpleTestUtils.createStructuralVariant(chr1, 1201, chr1, 1300, StructuralVariantType.DUP, 0.5, 0.5).build();
        mGermlineDeletions.findDeletions(List.of(pcn0, pcn1), List.of(or0, or1), List.of(svDup));
        List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        assertEquals(gd1_1.GeneName, deletions.get(0).GeneName);

        setup();
        driverMap = Map.of(gd1_1.GeneName, driver1_1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(chr1, List.of(gd1_1));
        mGermlineDeletions = new GermlineDeletions(driverMap, mEnsemblDataCache, mGermlineDeletionFrequency);

        // BND
        final StructuralVariant svBnd = PurpleTestUtils.createStructuralVariant(chr1, 1201, chr1, 1300, StructuralVariantType.BND, 0.5, 0.5).build();
        mGermlineDeletions.findDeletions(List.of(pcn0, pcn1), List.of(or0, or1), List.of(svBnd));
        deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        assertEquals(gd1_1.GeneName, deletions.get(0).GeneName);

        setup();
        driverMap = Map.of(gd1_1.GeneName, driver1_1);
        mEnsemblDataCache.transcriptData = List.of(td1_1);
        mEnsemblDataCache.chrGeneMap = Map.of(chr1, List.of(gd1_1));
        mGermlineDeletions = new GermlineDeletions(driverMap, mEnsemblDataCache, mGermlineDeletionFrequency);

        // INS
        final StructuralVariant svIns = PurpleTestUtils.createStructuralVariant(chr1, 1201, chr1, 1300, StructuralVariantType.INS, 0.5, 0.5).build();
        mGermlineDeletions.findDeletions(List.of(pcn0, pcn1), List.of(or0, or1), List.of(svIns));
        deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        assertEquals(gd1_1.GeneName, deletions.get(0).GeneName);
    }

    @Test
    public void multipleStructuralVariantsMatchingRegion()
    {
        final PurpleCopyNumber pcn0 = pcn(chr1, 1, 1000, 2);
        final PurpleCopyNumber pcn1 = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or0 = or(chr1, 1, 1000, GermlineStatus.DIPLOID, 1, 1, 1, 0.5);
        final ObservedRegion or1 = or(chr1, 1001, 2000, GermlineStatus.HOM_DELETION, 1001, 1001, 1, 0.5);

        // Create multiple SVs that match the region
        final StructuralVariant svDel = PurpleTestUtils.createStructuralVariant(chr1, 1201, chr1, 1300, StructuralVariantType.DEL, 0.5, 0.5).build();
        final StructuralVariant svInv = PurpleTestUtils.createStructuralVariant(chr1, 1250, chr1, 1350, StructuralVariantType.INV, 0.5, 0.5).build();
        final StructuralVariant svDup = PurpleTestUtils.createStructuralVariant(chr1, 1150, chr1, 1400, StructuralVariantType.DUP, 0.5, 0.5).build();

        // DEL type should be preferred over other types
        mGermlineDeletions.findDeletions(List.of(pcn0, pcn1), List.of(or0, or1), List.of(svInv, svDel, svDup));
        List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(1, deletions.size());
        assertEquals(svDel.start().position(), deletions.get(0).RegionStart);
        assertEquals(svDel.end().position(), deletions.get(0).RegionEnd);
    }

    @Test
    public void structuralVariantMatchingBothStartAndEnd()
    {
        final PurpleCopyNumber pcn0 = pcn(chr1, 1, 1000, 2);
        final PurpleCopyNumber pcn1 = pcn(chr1, 1001, 2000, 2);
        final PurpleCopyNumber pcn2 = pcn(chr1, 2001, 3000, 2);
        final ObservedRegion or0 = or(chr1, 1, 1000, GermlineStatus.DIPLOID, 1, 1, 1, 0.5);
        final ObservedRegion or1 = or(chr1, 1001, 2000, GermlineStatus.HOM_DELETION, 1001, 1001, 1, 0.5);
        final ObservedRegion or2 = or(chr1, 2001, 3000, GermlineStatus.DIPLOID, 2001, 2001, 1, 0.5);

        // Create an SV that matches both start and end of the region
        final StructuralVariant sv = PurpleTestUtils.createStructuralVariant(chr1, 1001, chr1, 2000, StructuralVariantType.DEL, 0.5, 0.5).build();

        mGermlineDeletions.findDeletions(List.of(pcn0, pcn1, pcn2), List.of(or0, or1, or2), List.of(sv));
        List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(2, deletions.size());

        // Both deletions should have the same start and end positions
        for (GermlineDeletion deletion : deletions) {
            assertEquals(sv.start().position(), deletion.RegionStart);
            assertEquals(sv.end().position(), deletion.RegionEnd);
        }

        // One deletion should be for each gene
        boolean foundGene1 = false;
        boolean foundGene2 = false;
        for (GermlineDeletion deletion : deletions) {
            if (deletion.GeneName.equals(gd1_1.GeneName)) {
                foundGene1 = true;
            } else if (deletion.GeneName.equals(gd1_2.GeneName)) {
                foundGene2 = true;
            }
        }

        assertTrue(foundGene1);
        assertTrue(foundGene2);
    }

    @Test
    public void hetDeletionGermlineStatus()
    {
        final PurpleCopyNumber pcn = pcn(chr1, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1000, 2001, GermlineStatus.HET_DELETION, 1000, 1000, 1, 0.5);
        mGermlineDeletions.findDeletions(List.of(pcn), List.of(or1), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(2, deletions.size());
        assertEquals(CommonVcfTags.PASS, deletions.get(0).Filter);
        assertEquals(ReportedStatus.REPORTED, deletions.get(0).Reported);
        assertEquals(GermlineStatus.HET_DELETION, deletions.get(0).NormalStatus);
    }

    @Test
    public void emptyCopyNumbersList()
    {
        final ObservedRegion or1 = or(chr1, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 1, 0.5);
        mGermlineDeletions.findDeletions(List.of(), List.of(or1), List.of());
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(2, deletions.size());
        assertEquals(CommonVcfTags.PASS, deletions.get(0).Filter);
        assertEquals(ReportedStatus.REPORTED, deletions.get(0).Reported);
    }

    @Test
    public void multipleChromosomesWithStructuralVariants()
    {
        final PurpleCopyNumber pcn1 = pcn(chr1, 1001, 2000, 2);
        final PurpleCopyNumber pcn2 = pcn(chr2, 1001, 2000, 2);
        final ObservedRegion or1 = or(chr1, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 1, 0.5);
        final ObservedRegion or2 = or(chr2, 1000, 2001, GermlineStatus.HOM_DELETION, 1000, 1000, 1, 0.5);

        final StructuralVariant sv1 = PurpleTestUtils.createStructuralVariant(chr1, 1201, chr1, 1300, StructuralVariantType.DEL, 0.5, 0.5).build();
        final StructuralVariant sv2 = PurpleTestUtils.createStructuralVariant(chr2, 1201, chr2, 1300, StructuralVariantType.DEL, 0.5, 0.5).build();

        mGermlineDeletions.findDeletions(List.of(pcn1, pcn2), List.of(or1, or2), List.of(sv1, sv2));
        final List<GermlineDeletion> deletions = mGermlineDeletions.getDeletions();
        assertEquals(2, deletions.size());

        // Check that we have one deletion from each chromosome
        boolean foundChr1 = false;
        boolean foundChr2 = false;
        for (GermlineDeletion deletion : deletions) {
            if (deletion.Chromosome.equals(chr1)) {
                foundChr1 = true;
                assertEquals(sv1.start().position(), deletion.RegionStart);
                assertEquals(sv1.end().position(), deletion.RegionEnd);
            } else if (deletion.Chromosome.equals(chr2)) {
                foundChr2 = true;
                assertEquals(sv2.start().position(), deletion.RegionStart);
                assertEquals(sv2.end().position(), deletion.RegionEnd);
            }
        }

        assertTrue(foundChr1);
        assertTrue(foundChr2);
    }

    private PurpleCopyNumber pcn(String chromosome, int start, int end, double majorAlleleCopyNumber)
    {
        final ImmutablePurpleCopyNumber result =
                PurpleTestUtils.createCopyNumber(chromosome, start, end, 2 * majorAlleleCopyNumber).build();
        Preconditions.checkState(result.majorAlleleCopyNumber() == majorAlleleCopyNumber);
        return result;
    }

    private ObservedRegion or(String chromosome,
            int start,
            int end,
            GermlineStatus status)
    {
        return or(chromosome, start, end, status, start, start, 2.0, 0.9);
    }

    private ObservedRegion or(String chromosome,
            int start,
            int end,
            GermlineStatus status,
            int minStart,
            int maxStart,
            double refNormalisedCopyNumber,
            double observedNormalRatio)
    {
        return new ObservedRegion(
                chromosome,
                start,
                end,
                false,
                null, // SegmentSupport
                0,   // BAF count
                0.0, // observed BAF
                0,   // depth window count
                0.0, // observed tumor ratio
                observedNormalRatio,
                0.0, // unnormalisedObservedNormalRatio
                status,
                false, // svCluster
                0.50, // gc content
                minStart,
                maxStart,
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
        TranscriptData result = new TranscriptData(transcriptId, "TD" + transcriptId, geneId, isCanonical, (byte) 0, 0,0,0,0,null,null);
        Preconditions.checkState(Objects.equals(result.GeneId, geneId));
        Preconditions.checkState(result.TransId==transcriptId);
        Preconditions.checkState(result.IsCanonical==isCanonical);
        List<ExonData> exonDataList = List.of(exons);
        exonDataList.forEach(exonData -> Preconditions.checkState(Objects.equals(exonData.TransId, transcriptId)));
        result.setExons(exonDataList);
        return result;
    }

    private DriverGene driverGene(String geneName, DriverGeneGermlineReporting reporting)
    {
        final ImmutableDriverGene result = ImmutableDriverGene.builder()
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
                .amplificationRatio(0.23)
                .reportSomaticHotspot(true)
                .likelihoodType(DriverCategory.TSG)
                .reportGermlineDisruption(DriverGeneGermlineReporting.ANY)
                .reportPGX(true)
                .build();
        Preconditions.checkState(Objects.equals(result.gene(), geneName));
        Preconditions.checkState(Objects.equals(result.reportGermlineDeletion(), reporting));
        return result;
    }
}
