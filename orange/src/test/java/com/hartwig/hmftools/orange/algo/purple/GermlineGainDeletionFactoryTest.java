package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.orange.TestDataUtils;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class GermlineGainDeletionFactoryTest
{
    private static final String TEST_GENE = "gene";
    private static final double EPSILON = 1.0E-2;

    @Test
    public void canTransformReportableHomDeletionsToPartial()
    {
        GermlineGainDeletionFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        GermlineAmpDel reportablePartialHom1 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 400, 700);
        GermlineAmpDel reportablePartialHom2 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 800, 1000);
        List<GermlineAmpDel> deletions = List.of(reportablePartialHom1, reportablePartialHom2);

        GeneCopyNumber partialLoss = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 1D, 4D);

        PurpleDriver purpleDriver = TestPurpleGainDeletionFactory.driverBuilder()
                .gene(TEST_GENE)
                .type(PurpleDriverType.GERMLINE_DELETION)
                .build();

        List<PurpleGainDeletion> gainDels = factory.createGermlineGainDeletions(deletions, List.of(purpleDriver), List.of(partialLoss));
        PurpleGainDeletion gainDel = gainDels.get(0);

        assertEquals(1, gainDels.size());
        assertEquals(ReportedStatus.REPORTED, gainDel.driver().reportedStatus());
        assertEquals(CopyNumberInterpretation.PARTIAL_DEL, gainDel.interpretation());
        assertEquals(TEST_GENE, gainDel.gene());
        assertEquals(0, gainDel.minCopies(), EPSILON);
        assertEquals(4, gainDel.maxCopies(), EPSILON);
    }

    @Test
    public void canTransformReportableHomDeletionToFull()
    {
        GermlineGainDeletionFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        GermlineAmpDel reportableFullHom =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 100, 1200);
        GeneCopyNumber fullLoss = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 1D, 1D);

        PurpleDriver purpleDriver = TestPurpleGainDeletionFactory.driverBuilder()
                .gene(TEST_GENE)
                .type(PurpleDriverType.GERMLINE_DELETION)
                .build();

        List<PurpleGainDeletion> gainDels = factory.createGermlineGainDeletions(List.of(reportableFullHom), List.of(purpleDriver), List.of(fullLoss));
        PurpleGainDeletion gainDel = gainDels.get(0);

        assertEquals(1, gainDels.size());
        assertEquals(ReportedStatus.REPORTED, gainDel.driver().reportedStatus());
        assertEquals(CopyNumberInterpretation.FULL_DEL, gainDel.interpretation());
        assertEquals(TEST_GENE, gainDel.gene());
        assertEquals(0, gainDel.minCopies(), EPSILON);
        assertEquals(0, gainDel.maxCopies(), EPSILON);
    }

    @Test
    public void canTransformReportablePartialHomDeletionsToFullGeneLoss()
    {
        GermlineGainDeletionFactory factory = createTestFactory();

        // Gene runs from 150 to 950
        // Exons are 250-350, 450-550 and 600-900
        GermlineAmpDel reportablePartial1 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.1D, 200, 300);
        GermlineAmpDel reportablePartial2 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 300, 500);
        GermlineAmpDel reportablePartial3 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0D, 500, 800);
        GermlineAmpDel reportablePartial4 =
                GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 0.2D, 700, 2000);
        List<GermlineAmpDel> deletions = List.of(reportablePartial1, reportablePartial2, reportablePartial3, reportablePartial4);

        GeneCopyNumber partialLoss = GeneCopyNumberTestFactory.createGeneCopyNumber(TEST_GENE, 1D, 4D);

        PurpleDriver purpleDriver = TestPurpleGainDeletionFactory.driverBuilder()
                .gene(TEST_GENE)
                .type(PurpleDriverType.GERMLINE_DELETION)
                .build();

        List<PurpleGainDeletion> gainDels = factory.createGermlineGainDeletions(deletions, List.of(purpleDriver), List.of(partialLoss));
        PurpleGainDeletion gainDel = gainDels.get(0);

        assertEquals(1, gainDels.size());
        assertEquals(ReportedStatus.REPORTED, gainDel.driver().reportedStatus());
        assertEquals(CopyNumberInterpretation.FULL_DEL, gainDel.interpretation());
        assertEquals(TEST_GENE, gainDel.gene());
        assertEquals(0, gainDel.minCopies(), EPSILON);
        assertEquals(0.2, gainDel.maxCopies(), EPSILON);
    }

    @NotNull
    private static GermlineGainDeletionFactory createTestFactory()
    {
        EnsemblDataCache ensemblDataCache = TestDataUtils.loadTestCache();
        return new GermlineGainDeletionFactory(ensemblDataCache);
    }
}