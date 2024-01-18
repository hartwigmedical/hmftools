package com.hartwig.hmftools.orange.algo.purple;

import static org.apache.commons.math3.util.Precision.EPSILON;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

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
    public void canCreateGermlineFullLosses()
    {
        // Gene is needed to be able to match with ensembl test data
        GermlineDeletion hetReported = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1);
        GermlineDeletion homReported = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 2);
        GermlineDeletion homUnreported = GermlineDeletionTestFactory.create(TEST_GENE, false, GermlineStatus.HOM_DELETION, 3);

        PurpleData purple = ImmutablePurpleData.builder()
                .from(PurpleTestFactory.createMinimalTestPurpleData())
                .addAllSomaticGeneCopyNumbers(GeneCopyNumberTestFactory.builder().chromosome("1").geneName(TEST_GENE).build())
                .allGermlineStructuralVariants(Lists.newArrayList())
                .addAllGermlineDeletions(hetReported, homReported, homUnreported)
                .addReportableGermlineDeletions(hetReported, homReported)
                .build();

        PurpleInterpreter interpreter = createRealInterpreter();
        PurpleRecord interpreted = interpreter.interpret(purple);
        assertEquals(1, interpreted.reportableGermlineFullLosses().size());
    }

    @Test
    public void canImplyDeletionsFromBreakendsWildtypeLostHomozygous()
    {
        LinxBreakend left = LinxOrangeTestFactory.breakendBuilder()
                .reportedDisruption(true)
                .gene(TEST_GENE)
                .transcriptId("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(0.3)
                .build();
        LinxBreakend right = LinxOrangeTestFactory.breakendBuilder()
                .reportedDisruption(true)
                .gene(TEST_GENE)
                .transcriptId("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(0.4)
                .build();
        DriverGene driverGene = createDriverGene(DriverGeneGermlineReporting.WILDTYPE_LOST);

        StructuralVariant shortSv = create("vcf id 1", 10, 20);
        LinxSvAnnotation svAnnotation = LinxOrangeTestFactory.svAnnotationBuilder().svId(1).vcfId(shortSv.id()).build();

        List<GermlineDeletion> impliedMatch = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(),
                Lists.newArrayList(left, right),
                Lists.newArrayList(shortSv),
                Lists.newArrayList(svAnnotation),
                Lists.newArrayList(driverGene));

        assertEquals(1, impliedMatch.size());
        GermlineDeletion deletion = impliedMatch.get(0);
        assertEquals(TEST_GENE, deletion.GeneName);
        assertEquals(GermlineStatus.HOM_DELETION, deletion.TumorStatus);
        assertEquals(0.4, deletion.TumorCopyNumber, EPSILON);
        assertEquals(10, deletion.RegionStart);
        assertEquals(20, deletion.RegionEnd);

        GermlineDeletion existingDeletion = GermlineDeletionTestFactory.create(TEST_GENE);
        List<GermlineDeletion> impliedExisting = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(existingDeletion),
                Lists.newArrayList(left, right),
                Lists.newArrayList(shortSv),
                Lists.newArrayList(svAnnotation),
                Lists.newArrayList(driverGene));

        assertEquals(0, impliedExisting.size());

        StructuralVariant longSv = create("vcf id 1", 10, 200000);
        List<GermlineDeletion> impliedTooLong = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(existingDeletion),
                Lists.newArrayList(left, right),
                Lists.newArrayList(longSv),
                Lists.newArrayList(svAnnotation),
                Lists.newArrayList(driverGene));

        assertEquals(0, impliedTooLong.size());
    }

    @Test
    public void canImplyDeletionsFromBreakendsWildtypeLostHeterozygous()
    {
        LinxBreakend left = LinxOrangeTestFactory.breakendBuilder()
                .reportedDisruption(true)
                .gene(TEST_GENE)
                .transcriptId("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(1.2)
                .build();
        LinxBreakend right = LinxOrangeTestFactory.breakendBuilder()
                .reportedDisruption(true)
                .gene(TEST_GENE)
                .transcriptId("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(0.9)
                .build();
        DriverGene driverGene = createDriverGene(DriverGeneGermlineReporting.WILDTYPE_LOST);

        StructuralVariant shortSv = create("vcf id 1", 10, 20);
        LinxSvAnnotation svAnnotation = LinxOrangeTestFactory.svAnnotationBuilder().svId(1).vcfId(shortSv.id()).build();

        List<GermlineDeletion> impliedMatch = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(),
                Lists.newArrayList(left, right),
                Lists.newArrayList(shortSv),
                Lists.newArrayList(svAnnotation),
                Lists.newArrayList(driverGene));

        assertEquals(0, impliedMatch.size());
    }

    @Test
    public void canImplyDeletionsFromBreakendsAnyHomozygous()
    {
        LinxBreakend left = LinxOrangeTestFactory.breakendBuilder()
                .reportedDisruption(true)
                .gene(TEST_GENE)
                .transcriptId("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(0.3)
                .build();
        LinxBreakend right = LinxOrangeTestFactory.breakendBuilder()
                .reportedDisruption(true)
                .gene(TEST_GENE)
                .transcriptId("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(0.4)
                .build();
        DriverGene driverGene = createDriverGene(DriverGeneGermlineReporting.ANY);

        StructuralVariant shortSv = create("vcf id 1", 10, 20);
        LinxSvAnnotation svAnnotation = LinxOrangeTestFactory.svAnnotationBuilder().svId(1).vcfId(shortSv.id()).build();

        List<GermlineDeletion> impliedMatch = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(),
                Lists.newArrayList(left, right),
                Lists.newArrayList(shortSv),
                Lists.newArrayList(svAnnotation),
                Lists.newArrayList(driverGene));

        assertEquals(1, impliedMatch.size());
        GermlineDeletion deletion = impliedMatch.get(0);
        assertEquals(TEST_GENE, deletion.GeneName);
        assertEquals(GermlineStatus.HOM_DELETION, deletion.TumorStatus);
        assertEquals(0.4, deletion.TumorCopyNumber, EPSILON);
        assertEquals(10, deletion.RegionStart);
        assertEquals(20, deletion.RegionEnd);

        GermlineDeletion existingDeletion = GermlineDeletionTestFactory.create(TEST_GENE);
        List<GermlineDeletion> impliedExisting = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(existingDeletion),
                Lists.newArrayList(left, right),
                Lists.newArrayList(shortSv),
                Lists.newArrayList(svAnnotation),
                Lists.newArrayList(driverGene));

        assertEquals(0, impliedExisting.size());

        StructuralVariant longSv = create("vcf id 1", 10, 200000);
        List<GermlineDeletion> impliedTooLong = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(existingDeletion),
                Lists.newArrayList(left, right),
                Lists.newArrayList(longSv),
                Lists.newArrayList(svAnnotation),
                Lists.newArrayList(driverGene));

        assertEquals(0, impliedTooLong.size());
    }

    @Test
    public void canImplyDeletionsFromBreakendsAnyHeterozygous()
    {
        LinxBreakend left = LinxOrangeTestFactory.breakendBuilder()
                .reportedDisruption(true)
                .gene(TEST_GENE)
                .transcriptId("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(1.2)
                .build();
        LinxBreakend right = LinxOrangeTestFactory.breakendBuilder()
                .reportedDisruption(true)
                .gene(TEST_GENE)
                .transcriptId("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(0.9)
                .build();
        DriverGene driverGene = createDriverGene(DriverGeneGermlineReporting.ANY);

        StructuralVariant shortSv = create("vcf id 1", 10, 20);
        LinxSvAnnotation svAnnotation = LinxOrangeTestFactory.svAnnotationBuilder().svId(1).vcfId(shortSv.id()).build();

        List<GermlineDeletion> impliedMatch = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(),
                Lists.newArrayList(left, right),
                Lists.newArrayList(shortSv),
                Lists.newArrayList(svAnnotation),
                Lists.newArrayList(driverGene));

        assertEquals(1, impliedMatch.size());
        GermlineDeletion deletion = impliedMatch.get(0);
        assertEquals(TEST_GENE, deletion.GeneName);
        assertEquals(GermlineStatus.HET_DELETION, deletion.TumorStatus);
        assertEquals(1.2, deletion.TumorCopyNumber, EPSILON);
        assertEquals(10, deletion.RegionStart);
        assertEquals(20, deletion.RegionEnd);

        GermlineDeletion existingDeletion = GermlineDeletionTestFactory.create(TEST_GENE);
        List<GermlineDeletion> impliedExisting = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(existingDeletion),
                Lists.newArrayList(left, right),
                Lists.newArrayList(shortSv),
                Lists.newArrayList(svAnnotation),
                Lists.newArrayList(driverGene));

        assertEquals(0, impliedExisting.size());

        StructuralVariant longSv = create("vcf id 1", 10, 200000);
        List<GermlineDeletion> impliedTooLong = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(existingDeletion),
                Lists.newArrayList(left, right),
                Lists.newArrayList(longSv),
                Lists.newArrayList(svAnnotation),
                Lists.newArrayList(driverGene));

        assertEquals(0, impliedTooLong.size());
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
        PurpleVariantFactory purpleVariantFactory = new PurpleVariantFactory(new PaveAlgo(ensemblDataCache));
        GermlineGainLossFactory germlineGainLossFactory = new GermlineGainLossFactory(ensemblDataCache);
        GermlineHeterozygousDeletionFactory germlineHeterozygousDeletionFactory = new GermlineHeterozygousDeletionFactory(ensemblDataCache);

        return new PurpleInterpreter(purpleVariantFactory,
                germlineGainLossFactory,
                germlineHeterozygousDeletionFactory,
                Lists.newArrayList(),
                TestLinxInterpretationFactory.createMinimalTestLinxData(),
                null);
    }

    @NotNull
    private static StructuralVariant create(@NotNull String vcfId, int start, int end)
    {
        return ImmutableStructuralVariantImpl.builder()
                .id(vcfId)
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
}