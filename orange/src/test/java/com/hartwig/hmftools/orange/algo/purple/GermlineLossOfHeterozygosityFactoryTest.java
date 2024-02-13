package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.datamodel.purple.GeneProportion;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineLossOfHeterozygosityFactoryTest
{
    private static final String TEST_GENE = "gene";
    private static final double EPSILON = 1.0E-2;

    @Test
    public void canFilterHomDeletion()
    {
        GermlineLossOfHeterozygosityFactory factory = createTestFactory();

        GermlineDeletion reportableHom = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION);
        assertTrue(factory.getReportabilityMap(Lists.newArrayList(reportableHom), Lists.newArrayList()).isEmpty());
    }

    @Test
    public void canTransformReportablePartialHetDeletionsToLOH()
    {
        GermlineLossOfHeterozygosityFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        GermlineDeletion reportablePartial1 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1D, 400, 700);
        GermlineDeletion reportablePartial2 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 2D, 800, 1000);
        List<GermlineDeletion> deletions = Lists.newArrayList(reportablePartial2, reportablePartial1);

        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).minCopyNumber(2D).maxCopyNumber(4D).build();

        Map<PurpleLossOfHeterozygosity, Boolean> map = factory.getReportabilityMap(deletions, Lists.newArrayList(geneCopyNumber));
        PurpleLossOfHeterozygosity loh = map.keySet().iterator().next();

        assertEquals(1, map.keySet().size());
        assertTrue(map.get(loh));
        assertEquals(GeneProportion.PARTIAL_GENE, loh.geneProportion());
        assertEquals(TEST_GENE, loh.gene());
        assertEquals(1, loh.minCopies(), EPSILON);
        assertEquals(4, loh.maxCopies(), EPSILON);
    }

    @Test
    public void canTransformNonReportableFullHetDeletionToLOH()
    {
        GermlineLossOfHeterozygosityFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        GermlineDeletion reportablePartial =
                GermlineDeletionTestFactory.create(TEST_GENE, false, GermlineStatus.HET_DELETION, 1D, 100, 1000);
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).minCopyNumber(2D).maxCopyNumber(4D).build();

        Map<PurpleLossOfHeterozygosity, Boolean> map =
                factory.getReportabilityMap(Lists.newArrayList(reportablePartial), Lists.newArrayList(geneCopyNumber));
        PurpleLossOfHeterozygosity loh = map.keySet().iterator().next();

        assertEquals(1, map.keySet().size());
        assertFalse(map.get(loh));
        assertEquals(GeneProportion.FULL_GENE, loh.geneProportion());
        assertEquals(TEST_GENE, loh.gene());
        assertEquals(1, loh.minCopies(), EPSILON);
        assertEquals(1, loh.maxCopies(), EPSILON);
    }

    @Test
    public void canTransformReportablePartialHetDeletionsToFullGeneLOH()
    {
        GermlineLossOfHeterozygosityFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        GermlineDeletion reportablePartial1 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 0.9D, 200, 375);
        GermlineDeletion reportablePartial2 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1.3D, 425, 450);
        GermlineDeletion reportablePartial3 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1.3D, 425, 525);
        GermlineDeletion reportablePartial4 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 3D, 500, 2000);
        List<GermlineDeletion> germlineDeletions =
                Lists.newArrayList(reportablePartial4, reportablePartial3, reportablePartial2, reportablePartial1);

        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).minCopyNumber(2D).maxCopyNumber(4D).build();

        Map<PurpleLossOfHeterozygosity, Boolean> map = factory.getReportabilityMap(germlineDeletions, Lists.newArrayList(geneCopyNumber));
        PurpleLossOfHeterozygosity loh = map.keySet().iterator().next();

        assertEquals(1, map.keySet().size());
        assertTrue(map.get(loh));
        assertEquals(GeneProportion.FULL_GENE, loh.geneProportion());
        assertEquals(TEST_GENE, loh.gene());
        assertEquals(0.9, loh.minCopies(), EPSILON);
        assertEquals(3, loh.maxCopies(), EPSILON);
    }

    @NotNull
    private static GermlineLossOfHeterozygosityFactory createTestFactory()
    {
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.loadTestCache();
        return new GermlineLossOfHeterozygosityFactory(ensemblDataCache);
    }
}