package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.purple.GeneCopyNumberTestFactory;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.GermlineDeletionTestFactory;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantImpl;
import com.hartwig.hmftools.common.sv.ImmutableStructuralVariantLegImpl;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.orange.algo.linx.LinxOrangeTestFactory;
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
    public void canCreateReportableGermlineFullDels()
    {
        // Gene is needed to be able to match with ensembl test data
        GermlineAmpDel hetReported = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HET_DELETION, 1);
        GermlineAmpDel homReported = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 2);

        PurpleData purple = createPurpleTestData(Lists.newArrayList(hetReported, homReported));

        PurpleInterpreter interpreter = createRealInterpreter();
        PurpleRecord interpreted = interpreter.interpret(purple);
        assertEquals(1, interpreted.germlineGainsDels().size());
    }

    @Test
    public void canHandleHalfReportableGermlineFullDels()
    {
        // Gene is needed to be able to match with ensembl test data
        GermlineAmpDel hetUnreported = GermlineDeletionTestFactory.create(TEST_GENE, false, GermlineStatus.HET_DELETION, 1);
        GermlineAmpDel homReported = GermlineDeletionTestFactory.create(TEST_GENE, true, GermlineStatus.HOM_DELETION, 2);

        PurpleData purple = createPurpleTestData(Lists.newArrayList(hetUnreported, homReported));

        PurpleInterpreter interpreter = createRealInterpreter();
        PurpleRecord interpreted = interpreter.interpret(purple);
        // assertEquals(1, interpreted.otherGermlineDeletions().size());
        assertEquals(1, interpreted.germlineGainsDels().size());
    }

    @Test
    public void canCreateNonReportableGermlineFullDels()
    {
        // Gene is needed to be able to match with ensembl test data
        GermlineAmpDel hetUnreported = GermlineDeletionTestFactory.create(TEST_GENE, false, GermlineStatus.HET_DELETION, 1);
        GermlineAmpDel homUnreported = GermlineDeletionTestFactory.create(TEST_GENE, false, GermlineStatus.HOM_DELETION, 2);

        PurpleData purple = createPurpleTestData(Lists.newArrayList(hetUnreported, homUnreported));

        PurpleInterpreter interpreter = createRealInterpreter();
        PurpleRecord interpreted = interpreter.interpret(purple);
        // assertEquals(1, interpreted.otherGermlineDeletions(). size());
        assertEquals(0, interpreted.germlineGainsDels().size());
    }


    @NotNull
    private static ImmutablePurpleData createPurpleTestData(@NotNull List<GermlineAmpDel> allGermlineDeletions)
    {
        return ImmutablePurpleData.builder()
                .from(PurpleTestFactory.createMinimalTestPurpleData())
                .addSomaticGeneCopyNumbers(GeneCopyNumberTestFactory.createGeneCopyNumber("1", TEST_GENE, 0, 0))
                .addAllGermlineDeletions(allGermlineDeletions)
                .germlineDeletions(allGermlineDeletions.stream().filter(d -> d.Reported == ReportedStatus.REPORTED).collect(Collectors.toList()))
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
        GermlineGainDeletionFactory germlineGainDeletionFactory = new GermlineGainDeletionFactory(ensemblDataCache);
        GermlineLossOfHeterozygosityFactory germlineLossOfHeterozygosityFactory = new GermlineLossOfHeterozygosityFactory(ensemblDataCache);

        return new PurpleInterpreter(
                purpleVariantFactory,
                germlineGainDeletionFactory,
                germlineLossOfHeterozygosityFactory);
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
                .hotspot(false)
                .build();
    }

    @NotNull
    private static LinxSvAnnotation createSvAnnotation(@NotNull StructuralVariant sv)
    {
        return LinxOrangeTestFactory.svAnnotationBuilder().svId(1).vcfId(sv.id()).build();
    }
}