package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

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
        assertTrue(factory.mapDeletions(Lists.newArrayList(reportableHom), Lists.newArrayList()).isEmpty());
    }

    @Test
    public void canTransformPartialHetDeletionToLOH()
    {
        GermlineLossOfHeterozygosityFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        GermlineDeletion reportablePartial =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1D, 400, 500);
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).minCopyNumber(2D).maxCopyNumber(4D).build();

        Map<PurpleLossOfHeterozygosity, GermlineDeletion> map =
                factory.mapDeletions(Lists.newArrayList(reportablePartial), Lists.newArrayList(geneCopyNumber));
        PurpleLossOfHeterozygosity heterozygousDeletion = map.keySet().iterator().next();

        assertEquals(reportablePartial, map.get(heterozygousDeletion));
        assertEquals(GeneProportion.PARTIAL_GENE, heterozygousDeletion.geneProportion());
        assertEquals(TEST_GENE, heterozygousDeletion.gene());
        assertEquals(1, heterozygousDeletion.minCopies(), EPSILON);
        assertEquals(4, heterozygousDeletion.maxCopies(), EPSILON);
    }

    @Test
    public void canTransformFullHetDeletionToLOH()
    {
        GermlineLossOfHeterozygosityFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        GermlineDeletion reportablePartial =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1D, 100, 1000);
        GeneCopyNumber geneCopyNumber = GeneCopyNumberTestFactory.builder().geneName(TEST_GENE).minCopyNumber(2D).maxCopyNumber(4D).build();

        Map<PurpleLossOfHeterozygosity, GermlineDeletion> map =
                factory.mapDeletions(Lists.newArrayList(reportablePartial), Lists.newArrayList(geneCopyNumber));
        PurpleLossOfHeterozygosity heterozygousDeletion = map.keySet().iterator().next();

        assertEquals(reportablePartial, map.get(heterozygousDeletion));
        assertEquals(GeneProportion.FULL_GENE, heterozygousDeletion.geneProportion());
        assertEquals(TEST_GENE, heterozygousDeletion.gene());
        assertEquals(1, heterozygousDeletion.minCopies(), EPSILON);
        assertEquals(1, heterozygousDeletion.maxCopies(), EPSILON);
    }

    @NotNull
    private static GermlineLossOfHeterozygosityFactory createTestFactory()
    {
        EnsemblDataCache ensemblDataCache = TestEnsemblDataCacheFactory.loadTestCache();
        return new GermlineLossOfHeterozygosityFactory(ensemblDataCache);
    }
}