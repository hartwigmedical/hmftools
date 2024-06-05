package com.hartwig.hmftools.orange.algo.immuno;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.orange.algo.linx.LinxOrangeTestFactory;
import com.hartwig.hmftools.orange.algo.linx.TestLinxInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGainLossFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGeneCopyNumberFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ImmuneEscapeInterpreterTest
{
    @Test
    public void doesNotCrashOnMinimalData()
    {
        PurpleRecord purple = TestPurpleInterpretationFactory.createMinimalTestPurpleData();
        LinxRecord linx = TestLinxInterpretationFactory.createMinimalTestLinxData();

        assertNotNull(ImmuneEscapeInterpreter.interpret(purple, linx));
    }

    @Test
    public void canDetectFullSpectrumOnRealData()
    {
        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .addAllSomaticVariants(TestPurpleVariantFactory.builder()
                        .gene("HLA-C")
                        .canonicalImpact(TestPurpleVariantFactory.impactBuilder().codingEffect(PurpleCodingEffect.MISSENSE).build())
                        .biallelic(true)
                        .subclonalLikelihood(0D)
                        .build())
                .addAllSomaticGainsLosses(TestPurpleGainLossFactory.builder()
                        .gene("TABBP")
                        .interpretation(CopyNumberInterpretation.FULL_LOSS)
                        .isCanonical(false)
                        .build())
                .addAllSomaticGainsLosses(TestPurpleGainLossFactory.builder()
                        .gene("CD274")
                        .interpretation(CopyNumberInterpretation.PARTIAL_GAIN)
                        .isCanonical(true)
                        .build())
                .addAllSomaticGainsLosses(TestPurpleGainLossFactory.builder()
                        .gene("SETDB1")
                        .interpretation(CopyNumberInterpretation.FULL_GAIN)
                        .isCanonical(true)
                        .build())
                .build();

        LinxRecord linx = TestLinxInterpretationFactory.builder()
                .addSomaticHomozygousDisruptions(LinxOrangeTestFactory.homozygousDisruptionBuilder()
                        .gene("IFNGR2")
                        .isCanonical(true)
                        .build())
                .build();

        ImmuneEscapeRecord immuneEscape = ImmuneEscapeInterpreter.interpret(purple, linx);

        assertTrue(immuneEscape.hasHlaEscape());
        assertFalse(immuneEscape.hasAntigenPresentationPathwayEscape());
        assertTrue(immuneEscape.hasIFNGammaPathwayEscape());
        assertFalse(immuneEscape.hasPDL1OverexpressionEscape());
        assertFalse(immuneEscape.hasCD58InactivationEscape());
        assertTrue(immuneEscape.hasEpigeneticSETDB1Escape());
    }

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

    @Test
    public void canDetectPDL1OverexpressionEscape()
    {
        String matchingGene = ImmuneEscapeInterpreter.PD_L1_GENES.iterator().next();

        assertTrue(runWithPurple(withAmplification(matchingGene)).hasPDL1OverexpressionEscape());

        assertFalse(runWithPurple(withDeletion(matchingGene)).hasPDL1OverexpressionEscape());
        assertFalse(runWithPurple(withAmplification("random gene")).hasPDL1OverexpressionEscape());
    }

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

    @Test
    public void canDetectEpigeneticSETDB1Escape()
    {
        String matchingGene = ImmuneEscapeInterpreter.EPIGENETIC_SETDB1_GENES.iterator().next();

        assertTrue(runWithPurple(withAmplification(matchingGene)).hasEpigeneticSETDB1Escape());

        assertFalse(runWithPurple(withDeletion(matchingGene)).hasEpigeneticSETDB1Escape());
        assertFalse(runWithPurple(withAmplification("random gene")).hasEpigeneticSETDB1Escape());
    }

    @NotNull
    private static ImmuneEscapeRecord runWithPurple(@NotNull PurpleRecord purple)
    {
        return ImmuneEscapeInterpreter.interpret(purple, TestLinxInterpretationFactory.createMinimalTestLinxData());
    }

    @NotNull
    private static ImmuneEscapeRecord runWithLinx(@NotNull LinxRecord linx)
    {
        return ImmuneEscapeInterpreter.interpret(TestPurpleInterpretationFactory.createMinimalTestPurpleData(), linx);
    }

    @NotNull
    private static PurpleRecord withLOH(@NotNull String gene)
    {
        return TestPurpleInterpretationFactory.builder().addAllSomaticGeneCopyNumbers(
                        TestPurpleGeneCopyNumberFactory.builder()
                                .gene(gene).minMinorAlleleCopyNumber(0D).minCopyNumber(1D).build())
                .build();
    }

    @NotNull
    private static PurpleRecord withoutLOH(@NotNull String gene)
    {
        return TestPurpleInterpretationFactory.builder().addAllSomaticGeneCopyNumbers(
                        TestPurpleGeneCopyNumberFactory.builder()
                                .gene(gene).minMinorAlleleCopyNumber(1D).minCopyNumber(2D).build())
                .build();
    }

    @NotNull
    private PurpleRecord withClonalVariant(@NotNull String gene, @NotNull PurpleCodingEffect codingEffect, boolean biallelic)
    {
        return withVariant(gene, codingEffect, biallelic, 0D);
    }

    @NotNull
    private PurpleRecord withSubclonalVariant(@NotNull String gene, @NotNull PurpleCodingEffect codingEffect, boolean biallelic)
    {
        return withVariant(gene, codingEffect, biallelic, 1D);
    }

    @NotNull
    private static PurpleRecord withVariant(@NotNull String gene, @NotNull PurpleCodingEffect codingEffect,
            boolean biallelic, double subclonalLikelihood)
    {
        return TestPurpleInterpretationFactory.builder().addAllSomaticVariants(
                        TestPurpleVariantFactory.builder()
                                .gene(gene)
                                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().codingEffect(codingEffect).build())
                                .biallelic(biallelic)
                                .subclonalLikelihood(subclonalLikelihood)
                                .build())
                .build();
    }

    @NotNull
    private static PurpleRecord withDeletion(@NotNull String gene)
    {
        return TestPurpleInterpretationFactory.builder()
                .addAllSomaticGainsLosses(TestPurpleGainLossFactory.builder()
                        .gene(gene)
                        .isCanonical(true)
                        .interpretation(CopyNumberInterpretation.FULL_LOSS)
                        .build())
                .build();
    }

    @NotNull
    private static PurpleRecord withAmplification(@NotNull String gene)
    {
        return TestPurpleInterpretationFactory.builder()
                .addAllSomaticGainsLosses(TestPurpleGainLossFactory.builder()
                        .gene(gene)
                        .isCanonical(true)
                        .interpretation(CopyNumberInterpretation.FULL_GAIN)
                        .build())
                .build();
    }

    @NotNull
    private static LinxRecord withHomozygousDisruption(@NotNull String gene)
    {
        return TestLinxInterpretationFactory.builder()
                .addSomaticHomozygousDisruptions(LinxOrangeTestFactory.homozygousDisruptionBuilder()
                        .gene(gene)
                        .isCanonical(true)
                        .build())
                .build();
    }
}