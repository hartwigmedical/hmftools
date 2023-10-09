package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;

import com.google.common.collect.Lists;
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
    public void canImplyDeletionsFromBreakends()
    {
        LinxBreakend left = LinxOrangeTestFactory.breakendBuilder()
                .reported(true)
                .gene(TEST_GENE)
                .transcript("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(0.3)
                .build();
        LinxBreakend right = LinxOrangeTestFactory.breakendBuilder()
                .reported(true)
                .gene(TEST_GENE)
                .transcript("trans 1")
                .svId(1)
                .type(LinxBreakendType.DEL)
                .undisruptedCopyNumber(0.4)
                .build();

        StructuralVariant shortSv = create("vcf id 1", 10, 20);
        LinxSvAnnotation svAnnotation = LinxOrangeTestFactory.svAnnotationBuilder().svId(1).vcfId(shortSv.id()).build();

        List<GermlineDeletion> impliedMatch = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(),
                Lists.newArrayList(left, right),
                Lists.newArrayList(shortSv),
                Lists.newArrayList(svAnnotation));

        assertEquals(1, impliedMatch.size());
        GermlineDeletion deletion = impliedMatch.get(0);
        assertEquals(TEST_GENE, deletion.GeneName);
        assertEquals(GermlineStatus.HOM_DELETION, deletion.TumorStatus);

        GermlineDeletion existingDeletion = GermlineDeletionTestFactory.create(TEST_GENE);
        List<GermlineDeletion> impliedExisting = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(existingDeletion),
                Lists.newArrayList(left, right),
                Lists.newArrayList(shortSv),
                Lists.newArrayList(svAnnotation));

        assertEquals(0, impliedExisting.size());

        StructuralVariant longSv = create("vcf id 1", 10, 200000);
        List<GermlineDeletion> impliedTooLong = PurpleInterpreter.implyDeletionsFromBreakends(Lists.newArrayList(existingDeletion),
                Lists.newArrayList(left, right),
                Lists.newArrayList(longSv),
                Lists.newArrayList(svAnnotation));

        assertEquals(0, impliedTooLong.size());
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

        return new PurpleInterpreter(purpleVariantFactory,
                germlineGainLossFactory,
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