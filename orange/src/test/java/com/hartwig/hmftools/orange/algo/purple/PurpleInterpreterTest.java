package com.hartwig.hmftools.orange.algo.purple;

import static org.apache.commons.math3.util.Precision.EPSILON;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantImpl;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.orange.algo.linx.LinxOrangeTestFactory;
import com.hartwig.hmftools.orange.algo.linx.TestLinxInterpretationFactory;
import com.hartwig.hmftools.orange.algo.pave.PaveAlgo;
import com.hartwig.hmftools.orange.algo.pave.TestEnsemblDataCacheFactory;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PurpleInterpreterTest
{
    private static final String TEST_GENE = "gene";

    @Test
    public void canInterpretMinimalPurpleData()
    {
        PurpleInterpreter interpreter = createTestInterpreter();
        assertNotNull(interpreter.interpret(PurpleTestFactory.createMinimalTestPurpleData()));
    }

    @Test
    public void canCreateReportableGermlineFullLosses()
    {
        // Gene is needed to be able to match with ensembl test data
        GermlineDeletion hetReported = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1);
        GermlineDeletion homReported = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 2);

        PurpleData purple = createPurpleTestData(Lists.newArrayList(hetReported, homReported));

        PurpleInterpreter interpreter = createRealInterpreter();
        PurpleRecord interpreted = interpreter.interpret(purple);
        assertEquals(1, interpreted.allGermlineFullLosses().size());
        assertEquals(1, interpreted.reportableGermlineFullLosses().size());
        assertEquals(1, interpreted.allGermlineLossOfHeterozygosities().size());
        assertEquals(1, interpreted.reportableGermlineLossOfHeterozygosities().size());
    }

    @Test
    public void canHandleHalfReportableGermlineFullLosses()
    {
        // Gene is needed to be able to match with ensembl test data
        GermlineDeletion hetUnreported = GermlineDeletionTestFactory.create(TEST_GENE, false, GermlineStatus.HET_DELETION, 1);
        GermlineDeletion homReported = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 2);

        PurpleData purple = createPurpleTestData(Lists.newArrayList(hetUnreported, homReported));

        PurpleInterpreter interpreter = createRealInterpreter();
        PurpleRecord interpreted = interpreter.interpret(purple);
        assertEquals(1, interpreted.allGermlineFullLosses().size());
        assertEquals(1, interpreted.reportableGermlineFullLosses().size());
        assertEquals(1, interpreted.allGermlineLossOfHeterozygosities().size());
        assertEquals(0, interpreted.reportableGermlineLossOfHeterozygosities().size());
    }

    @Test
    public void canCreateNonReportableGermlineFullLosses()
    {
        // Gene is needed to be able to match with ensembl test data
        GermlineDeletion hetUnreported = GermlineDeletionTestFactory.create(TEST_GENE, false, GermlineStatus.HET_DELETION, 1);
        GermlineDeletion homUnreported = GermlineDeletionTestFactory.create(TEST_GENE, false, GermlineStatus.HOM_DELETION, 2);

        PurpleData purple = createPurpleTestData(Lists.newArrayList(hetUnreported, homUnreported));

        PurpleInterpreter interpreter = createRealInterpreter();
        PurpleRecord interpreted = interpreter.interpret(purple);
        assertEquals(1, interpreted.allGermlineFullLosses().size());
        assertEquals(0, interpreted.reportableGermlineFullLosses().size());
        assertEquals(1, interpreted.allGermlineLossOfHeterozygosities().size());
        assertEquals(0, interpreted.reportableGermlineLossOfHeterozygosities().size());
    }

    @Test
    public void canImplyDeletionFromBreakendsWhenHomozygousInTumorWithWildtypeLost()
    {
        double leftBreakendUndisruptedCopyNumber = 0.3;
        double rightBreakendUndisruptedCopyNumber = 0.4;
        DriverGeneGermlineReporting reportGermlineDeletion = DriverGeneGermlineReporting.WILDTYPE_LOST;
        int shortSvStartPosition = 10;
        int shortSvEndPosition = 20;

        List<GermlineDeletion> impliedMatch = doTestRunImplyDeletionsFromBreakends(
                leftBreakendUndisruptedCopyNumber,
                rightBreakendUndisruptedCopyNumber,
                shortSvStartPosition,
                shortSvEndPosition,
                reportGermlineDeletion,
                Lists.newArrayList()
        );

        assertEquals(1, impliedMatch.size());
        GermlineDeletion deletion = impliedMatch.get(0);
        assertEquals(TEST_GENE, deletion.GeneName);
        assertEquals(GermlineStatus.HOM_DELETION, deletion.TumorStatus);
        assertEquals(0.4, deletion.TumorCopyNumber, EPSILON);
        assertEquals(10, deletion.RegionStart);
        assertEquals(20, deletion.RegionEnd);

        GermlineDeletion existingDeletion = GermlineDeletionTestFactory.create(TEST_GENE);
        List<GermlineDeletion> impliedExisting = doTestRunImplyDeletionsFromBreakends(
                leftBreakendUndisruptedCopyNumber,
                rightBreakendUndisruptedCopyNumber,
                shortSvStartPosition,
                shortSvEndPosition,
                reportGermlineDeletion,
                Lists.newArrayList(existingDeletion)
        );

        assertEquals(0, impliedExisting.size());

        int longSvStartPosition = 10;
        int longSvEndPosition = 200000;
        List<GermlineDeletion> impliedTooLong = doTestRunImplyDeletionsFromBreakends(
                leftBreakendUndisruptedCopyNumber,
                rightBreakendUndisruptedCopyNumber,
                longSvStartPosition,
                longSvEndPosition,
                reportGermlineDeletion,
                Lists.newArrayList()
        );

        assertEquals(0, impliedTooLong.size());
    }

    @Test
    public void cannotImplyDeletionFromBreakendsWhenHeterozygousInTumorWithWildtypeLost()
    {
        double leftBreakendUndisruptedCopyNumber = 1.2;
        double rightBreakendUndisruptedCopyNumber = 0.9;
        DriverGeneGermlineReporting reportGermlineDeletion = DriverGeneGermlineReporting.WILDTYPE_LOST;
        int shortSvStartPosition = 10;
        int shortSvEndPosition = 20;

        List<GermlineDeletion> impliedMatch = doTestRunImplyDeletionsFromBreakends(
                leftBreakendUndisruptedCopyNumber,
                rightBreakendUndisruptedCopyNumber,
                shortSvStartPosition,
                shortSvEndPosition,
                reportGermlineDeletion,
                Lists.newArrayList()
        );

        assertEquals(0, impliedMatch.size());
    }

    @Test
    public void canImplyDeletionFromBreakendsWhenHomozygousInTumorWithAny()
    {
        double leftBreakendUndisruptedCopyNumber = 0.3;
        double rightBreakendUndisruptedCopyNumber = 0.4;
        DriverGeneGermlineReporting reportGermlineDeletion = DriverGeneGermlineReporting.ANY;
        int shortSvStartPosition = 10;
        int shortSvEndPosition = 20;

        List<GermlineDeletion> impliedMatch = doTestRunImplyDeletionsFromBreakends(
                leftBreakendUndisruptedCopyNumber,
                rightBreakendUndisruptedCopyNumber,
                shortSvStartPosition,
                shortSvEndPosition,
                reportGermlineDeletion,
                Lists.newArrayList()
        );

        assertEquals(1, impliedMatch.size());
        GermlineDeletion deletion = impliedMatch.get(0);
        assertEquals(TEST_GENE, deletion.GeneName);
        assertEquals(GermlineStatus.HOM_DELETION, deletion.TumorStatus);
        assertEquals(0.4, deletion.TumorCopyNumber, EPSILON);
        assertEquals(10, deletion.RegionStart);
        assertEquals(20, deletion.RegionEnd);

        GermlineDeletion existingDeletion = GermlineDeletionTestFactory.create(TEST_GENE);
        List<GermlineDeletion> impliedExisting = doTestRunImplyDeletionsFromBreakends(
                leftBreakendUndisruptedCopyNumber,
                rightBreakendUndisruptedCopyNumber,
                shortSvStartPosition,
                shortSvEndPosition,
                reportGermlineDeletion,
                Lists.newArrayList(existingDeletion)
        );

        assertEquals(0, impliedExisting.size());

        int longSvStartPosition = 10;
        int longSvEndPosition = 200000;
        List<GermlineDeletion> impliedTooLong = doTestRunImplyDeletionsFromBreakends(
                leftBreakendUndisruptedCopyNumber,
                rightBreakendUndisruptedCopyNumber,
                longSvStartPosition,
                longSvEndPosition,
                reportGermlineDeletion,
                Lists.newArrayList()
        );

        assertEquals(0, impliedTooLong.size());
    }

    @Test
    public void canImplyDeletionFromBreakendsWhenHeterozygousInTumorWithAny()
    {
        double leftBreakendUndisruptedCopyNumber = 1.2;
        double rightBreakendUndisruptedCopyNumber = 0.9;
        DriverGeneGermlineReporting reportGermlineDeletion = DriverGeneGermlineReporting.ANY;
        int shortSvStartPosition = 10;
        int shortSvEndPosition = 20;

        List<GermlineDeletion> impliedMatch = doTestRunImplyDeletionsFromBreakends(
                leftBreakendUndisruptedCopyNumber,
                rightBreakendUndisruptedCopyNumber,
                shortSvStartPosition,
                shortSvEndPosition,
                reportGermlineDeletion,
                Lists.newArrayList()
        );

        assertEquals(1, impliedMatch.size());
        GermlineDeletion deletion = impliedMatch.get(0);
        assertEquals(TEST_GENE, deletion.GeneName);
        assertEquals(GermlineStatus.HET_DELETION, deletion.TumorStatus);
        assertEquals(leftBreakendUndisruptedCopyNumber, deletion.TumorCopyNumber, EPSILON);
        assertEquals(shortSvStartPosition, deletion.RegionStart);
        assertEquals(shortSvEndPosition, deletion.RegionEnd);

        GermlineDeletion existingDeletion = GermlineDeletionTestFactory.create(TEST_GENE);
        List<GermlineDeletion> impliedExisting = doTestRunImplyDeletionsFromBreakends(
                leftBreakendUndisruptedCopyNumber,
                rightBreakendUndisruptedCopyNumber,
                shortSvStartPosition,
                shortSvEndPosition,
                reportGermlineDeletion,
                Lists.newArrayList(existingDeletion)
        );

        assertEquals(0, impliedExisting.size());

        int longSvStartPosition = 10;
        int longSvEndPosition = 200000;
        List<GermlineDeletion> impliedTooLong = doTestRunImplyDeletionsFromBreakends(
                leftBreakendUndisruptedCopyNumber,
                rightBreakendUndisruptedCopyNumber,
                longSvStartPosition,
                longSvEndPosition,
                reportGermlineDeletion,
                Lists.newArrayList()
        );

        assertEquals(0, impliedTooLong.size());
    }

    @NotNull
    private static ImmutablePurpleData createPurpleTestData(@NotNull List<GermlineDeletion> allGermlineDeletions)
    {
        return ImmutablePurpleData.builder()
                .from(PurpleTestFactory.createMinimalTestPurpleData())
                .addAllSomaticGeneCopyNumbers(GeneCopyNumberTestFactory.builder().chromosome("1").geneName(TEST_GENE).build())
                .allGermlineStructuralVariants(Lists.newArrayList())
                .allGermlineDeletions(allGermlineDeletions)
                .reportableGermlineDeletions(allGermlineDeletions.stream().filter(d -> d.Reported).collect(Collectors.toList()))
                .build();
    }

    @NotNull
    private static List<GermlineDeletion> doTestRunImplyDeletionsFromBreakends(double leftBreakendUndisruptedCopyNumber,
            double rightBreakendUndisruptedCopyNumber, int svStartPosition, int svEndPosition,
            DriverGeneGermlineReporting reportGermlineDeletion, List<GermlineDeletion> existingDeletions)
    {
        List<LinxBreakend> breakends = createBreakends(leftBreakendUndisruptedCopyNumber, rightBreakendUndisruptedCopyNumber);
        StructuralVariant sv = createStructuralVariant(svStartPosition, svEndPosition);
        List<DriverGene> driverGenes = Lists.newArrayList(createDriverGene(reportGermlineDeletion));

        return PurpleInterpreter.implyDeletionsFromBreakends(
                existingDeletions,
                breakends,
                Lists.newArrayList(sv),
                Lists.newArrayList(createSvAnnotation(sv)),
                driverGenes
        );
    }

    @NotNull
    private static List<LinxBreakend> createBreakends(double leftUndisruptedCopyNumber, double rightUndisruptedCopyNumber)
    {
        LinxBreakend left = LinxOrangeTestFactory.breakendBuilder()
                .reported(true)
                .gene(TEST_GENE)
                .transcript("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(leftUndisruptedCopyNumber)
                .build();
        LinxBreakend right = LinxOrangeTestFactory.breakendBuilder()
                .reported(true)
                .gene(TEST_GENE)
                .transcript("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(rightUndisruptedCopyNumber)
                .build();
        return Lists.newArrayList(left, right);
    }

    @NotNull
    private static ImmutableDriverGene createDriverGene(@NotNull DriverGeneGermlineReporting reportGermlineDeletion)
    {
        return ImmutableDriverGene.builder()
                .gene(TEST_GENE)
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(true)
                .reportDisruption(false)
                .reportAmplification(true)
                .reportSomaticHotspot(false)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .reportGermlineDisruption(DriverGeneGermlineReporting.ANY)
                .reportGermlineDeletion(reportGermlineDeletion)
                .likelihoodType(DriverCategory.ONCO)
                .reportPGX(false)
                .build();
    }

    @NotNull
    private static PurpleInterpreter createTestInterpreter()
    {
        return createInterpreter(TestEnsemblDataCacheFactory.createDummyCache());
    }

    @NotNull
    private static PurpleInterpreter createRealInterpreter()
    {
        return createInterpreter(TestEnsemblDataCacheFactory.loadTestCache());
    }

    @NotNull
    private static PurpleInterpreter createInterpreter(@NotNull EnsemblDataCache ensemblDataCache)
    {
        PaveAlgo pave = new PaveAlgo(ensemblDataCache, false);
        PurpleVariantFactory purpleVariantFactory = new PurpleVariantFactory(pave);
        GermlineGainLossFactory germlineGainLossFactory = new GermlineGainLossFactory(ensemblDataCache);
        GermlineLossOfHeterozygosityFactory germlineLossOfHeterozygosityFactory = new GermlineLossOfHeterozygosityFactory(ensemblDataCache);

        return new PurpleInterpreter(
                purpleVariantFactory,
                germlineGainLossFactory,
                germlineLossOfHeterozygosityFactory,
                Lists.newArrayList(),
                TestLinxInterpretationFactory.createMinimalTestLinxData(),
                null
        );
    }

    @NotNull
    private static StructuralVariant createStructuralVariant(int start, int end)
    {
        return ImmutableStructuralVariantImpl.builder()
                .id("vcf id 1")
                .start(ImmutableStructuralVariantLegImpl.builder()
                        .orientation((byte) 0)
                        .homology(Strings.EMPTY)
                        .anchoringSupportDistance(0)
                        .chromosome(Strings.EMPTY)
                        .position(start)
                        .build())
                .end(ImmutableStructuralVariantLegImpl.builder()
                        .orientation((byte) 0)
                        .homology(Strings.EMPTY)
                        .anchoringSupportDistance(0)
                        .chromosome(Strings.EMPTY)
                        .position(end)
                        .build())
                .insertSequence(Strings.EMPTY)
                .type(StructuralVariantType.DEL)
                .qualityScore(0D)
                .recovered(false)
                .hotspot(false)
                .build();
    }

    @NotNull
    private static LinxSvAnnotation createSvAnnotation(@NotNull StructuralVariant sv)
    {
        return LinxOrangeTestFactory.svAnnotationBuilder().svId(1).vcfId(sv.id()).build();
    }
}