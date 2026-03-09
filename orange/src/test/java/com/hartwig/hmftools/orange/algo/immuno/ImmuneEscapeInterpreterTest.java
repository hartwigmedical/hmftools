package com.hartwig.hmftools.orange.algo.immuno;

import static com.hartwig.hmftools.orange.algo.purple.PurpleTestFactory.purpleDriverBuilder;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.Collections;

import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.orange.algo.linx.TestLinxFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGainDeletionFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGeneCopyNumberFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;

import org.junit.Ignore;
import org.junit.Test;

public class ImmuneEscapeInterpreterTest
{
    @Test
    public void doesNotCrashOnMinimalData()
    {
        PurpleRecord purple = TestPurpleInterpretationFactory.createMinimalTestPurpleData();
        LinxRecord linx = TestLinxFactory.createMinimalTestLinxData();

        assertNotNull(ImmuneEscapeInterpreter.interpret(purple, linx));
    }

    @Ignore
    @Test
    public void canDetectFullSpectrumOnRealData()
    {
        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .addSomaticGainsDels(TestPurpleGainDeletionFactory.builder()
                        .driver(TestPurpleGainDeletionFactory.driverBuilder()
                                .gene("TABBP")
                                .isCanonical(false)
                                .build())
                        //.interpretation(CopyNumberInterpretation.FULL_DEL)
                        .build())
                .addSomaticGainsDels(TestPurpleGainDeletionFactory.builder()
                        .driver(TestPurpleGainDeletionFactory.driverBuilder()
                                .gene("CD274")
                                .isCanonical(true)
                                .build())
                        //.interpretation(CopyNumberInterpretation.PARTIAL_GAIN)
                        .build())
                .addSomaticGainsDels(TestPurpleGainDeletionFactory.builder()
                        .driver(TestPurpleGainDeletionFactory.driverBuilder()
                                .gene("SETDB1")
                                .isCanonical(true)
                                .build())
                        // .interpretation(CopyNumberInterpretation.FULL_GAIN)
                        .build())
                .build();

        LinxRecord linx = TestLinxFactory.linxRecordBuilder()
                .build();

        ImmuneEscapeRecord immuneEscape = ImmuneEscapeInterpreter.interpret(purple, linx);

        assertTrue(immuneEscape.hasHlaEscape());
        assertFalse(immuneEscape.hasAntigenPresentationPathwayEscape());
        assertTrue(immuneEscape.hasIFNGammaPathwayEscape());
        assertFalse(immuneEscape.hasPDL1OverexpressionEscape());
        assertFalse(immuneEscape.hasCD58InactivationEscape());
        assertTrue(immuneEscape.hasEpigeneticSETDB1Escape());
    }

    @Ignore
    @Test
    public void canDetectHlaEscape()
    {
        String matchingGene = ImmuneEscapeInterpreter.HLA_GENES.iterator().next();

        assertTrue(runWithPurple(withLOH(matchingGene)).hasHlaEscape());
        assertTrue(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false)).hasHlaEscape());
        assertTrue(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.MISSENSE, true)).hasHlaEscape());
        assertTrue(runWithPurple(withDeletion(matchingGene)).hasHlaEscape());
        assertTrue(runWithLinx(withHomozygousDisruption(matchingGene)).hasHlaEscape());

        assertFalse(runWithPurple(withoutLOH(matchingGene)).hasHlaEscape());
        assertFalse(runWithPurple(withSubclonalVariant(matchingGene, PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false)).hasHlaEscape());
        assertFalse(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.MISSENSE, false)).hasHlaEscape());
        assertFalse(runWithPurple(withAmplification(matchingGene)).hasHlaEscape());
        assertFalse(runWithPurple(withLOH("random gene")).hasHlaEscape());
    }

    @Ignore
    @Test
    public void canDetectAntigenPresentationPathwayEscape()
    {
        String matchingGene = ImmuneEscapeInterpreter.ANTIGEN_PRESENTATION_GENES.iterator().next();

        assertTrue(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false))
                .hasAntigenPresentationPathwayEscape());
        assertTrue(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.MISSENSE, true)).hasAntigenPresentationPathwayEscape());
        assertTrue(runWithPurple(withDeletion(matchingGene)).hasAntigenPresentationPathwayEscape());
        assertTrue(runWithLinx(withHomozygousDisruption(matchingGene)).hasAntigenPresentationPathwayEscape());

        assertFalse(runWithPurple(withLOH(matchingGene)).hasAntigenPresentationPathwayEscape());
        assertFalse(runWithPurple(withSubclonalVariant(matchingGene, PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false))
                .hasAntigenPresentationPathwayEscape());
        assertFalse(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.MISSENSE, false)).hasAntigenPresentationPathwayEscape());
        assertFalse(runWithPurple(withAmplification(matchingGene)).hasAntigenPresentationPathwayEscape());
        assertFalse(runWithPurple(withClonalVariant("random gene", PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false))
                .hasAntigenPresentationPathwayEscape());
    }

    @Ignore
    @Test
    public void canDetectIFNGammaPathwayEscape()
    {
        String matchingGene = ImmuneEscapeInterpreter.IFN_GAMMA_PATHWAY_GENES.iterator().next();

        assertTrue(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false))
                .hasIFNGammaPathwayEscape());
        assertTrue(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.MISSENSE, true)).hasIFNGammaPathwayEscape());
        assertTrue(runWithPurple(withDeletion(matchingGene)).hasIFNGammaPathwayEscape());
        assertTrue(runWithLinx(withHomozygousDisruption(matchingGene)).hasIFNGammaPathwayEscape());

        assertFalse(runWithPurple(withLOH(matchingGene)).hasIFNGammaPathwayEscape());
        assertFalse(runWithPurple(withSubclonalVariant(matchingGene, PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false))
                .hasIFNGammaPathwayEscape());
        assertFalse(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.MISSENSE, false)).hasIFNGammaPathwayEscape());
        assertFalse(runWithPurple(withAmplification(matchingGene)).hasIFNGammaPathwayEscape());
        assertFalse(runWithPurple(withClonalVariant("random gene", PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false))
                .hasIFNGammaPathwayEscape());
    }

    @Ignore
    @Test
    public void canDetectPDL1OverexpressionEscape()
    {
        String matchingGene = ImmuneEscapeInterpreter.PD_L1_GENES.iterator().next();

        assertTrue(runWithPurple(withAmplification(matchingGene)).hasPDL1OverexpressionEscape());

        assertFalse(runWithPurple(withDeletion(matchingGene)).hasPDL1OverexpressionEscape());
        assertFalse(runWithPurple(withAmplification("random gene")).hasPDL1OverexpressionEscape());
    }

    @Ignore
    @Test
    public void canDetectCD58InactivationEscape()
    {
        String matchingGene = ImmuneEscapeInterpreter.CD58_GENES.iterator().next();

        assertTrue(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false))
                .hasCD58InactivationEscape());
        assertTrue(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.MISSENSE, true)).hasCD58InactivationEscape());
        assertTrue(runWithPurple(withDeletion(matchingGene)).hasCD58InactivationEscape());
        assertTrue(runWithLinx(withHomozygousDisruption(matchingGene)).hasCD58InactivationEscape());

        assertFalse(runWithPurple(withLOH(matchingGene)).hasCD58InactivationEscape());
        assertFalse(runWithPurple(withSubclonalVariant(matchingGene, PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false))
                .hasCD58InactivationEscape());
        assertFalse(runWithPurple(withClonalVariant(matchingGene, PurpleCodingEffect.MISSENSE, false)).hasCD58InactivationEscape());
        assertFalse(runWithPurple(withAmplification(matchingGene)).hasCD58InactivationEscape());
        assertFalse(runWithPurple(withClonalVariant("random gene", PurpleCodingEffect.NONSENSE_OR_FRAMESHIFT, false))
                .hasCD58InactivationEscape());
    }

    @Ignore
    @Test
    public void canDetectEpigeneticSETDB1Escape()
    {
        String matchingGene = ImmuneEscapeInterpreter.EPIGENETIC_SETDB1_GENES.iterator().next();

        assertTrue(runWithPurple(withAmplification(matchingGene)).hasEpigeneticSETDB1Escape());

        assertFalse(runWithPurple(withDeletion(matchingGene)).hasEpigeneticSETDB1Escape());
        assertFalse(runWithPurple(withAmplification("random gene")).hasEpigeneticSETDB1Escape());
    }

    private static ImmuneEscapeRecord runWithPurple(final PurpleRecord purple)
    {
        return ImmuneEscapeInterpreter.interpret(purple, TestLinxFactory.createMinimalTestLinxData());
    }

    private static ImmuneEscapeRecord runWithLinx(final LinxRecord linx)
    {
        return ImmuneEscapeInterpreter.interpret(TestPurpleInterpretationFactory.createMinimalTestPurpleData(), linx);
    }

    private static PurpleRecord withLOH(final String gene)
    {
        return TestPurpleInterpretationFactory.builder().addSomaticDrivers(
                        purpleDriverBuilder().gene(gene).type(PurpleDriverType.LOH).build())
                .build();
    }

    private static PurpleRecord withoutLOH(final String gene)
    {
        return TestPurpleInterpretationFactory.builder().addSomaticDrivers(
                        purpleDriverBuilder().gene(gene).type(PurpleDriverType.AMP).build())
                .build();
    }

    private PurpleRecord withClonalVariant(final String gene, final PurpleCodingEffect codingEffect, boolean biallelic)
    {
        return withVariant(gene, codingEffect, biallelic, 0D);
    }

    private PurpleRecord withSubclonalVariant(final String gene, final PurpleCodingEffect codingEffect, boolean biallelic)
    {
        return withVariant(gene, codingEffect, biallelic, 1D);
    }

    private static PurpleRecord withVariant(final String gene, final PurpleCodingEffect codingEffect,
            boolean biallelic, double subclonalLikelihood)
    {
        return TestPurpleInterpretationFactory.builder().addSomaticVariants(
                        TestPurpleVariantFactory.builder()
                                .gene(gene)
                                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().codingEffect(codingEffect).build())
                                .biallelic(biallelic)
                                .subclonalLikelihood(subclonalLikelihood)
                                .build())
                .build();
    }

    private static PurpleRecord withDeletion(final String gene)
    {
        return TestPurpleInterpretationFactory.builder()
                .addSomaticGainsDels(TestPurpleGainDeletionFactory.builder()
                        .driver(TestPurpleGainDeletionFactory.driverBuilder()
                                .gene(gene)
                                .isCanonical(true)
                                .build())
                        // .interpretation(CopyNumberInterpretation.FULL_DEL)
                        .build())
                .build();
    }

    private static PurpleRecord withAmplification(final String gene)
    {
        return TestPurpleInterpretationFactory.builder()
                .addSomaticGainsDels(TestPurpleGainDeletionFactory.builder()
                        .driver(TestPurpleGainDeletionFactory.driverBuilder()
                                .gene(gene)
                                .isCanonical(true)
                                .build())
                        // .interpretation(CopyNumberInterpretation.FULL_GAIN)
                        .build())
                .build();
    }

    private static LinxRecord withHomozygousDisruption(final String gene)
    {
        return TestLinxFactory.linxRecordBuilder()
                .build();
    }
}