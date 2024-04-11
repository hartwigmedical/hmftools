package com.hartwig.hmftools.orange.algo.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.orange.TestOrangeReportFactory;
import com.hartwig.hmftools.orange.algo.linx.LinxOrangeTestFactory;
import com.hartwig.hmftools.orange.algo.linx.TestLinxInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGainLossFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGeneCopyNumberFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleLossOfHeterozygosityFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleQCFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class GermlineConversionTest
{
    private static final double EPSILON = 1.0E-10;
    private static final String TEST_GENE1 = "gene1";
    private static final String TEST_GENE2 = "gene2";
    private static final String TEST_TRANSCRIPT1 = "transcript1";
    private static final String TEST_TRANSCRIPT2 = "transcript2";

    @Test
    public void wipesAllGermlineUponConversion()
    {
        OrangeRecord report = TestOrangeReportFactory.builder()
                .germlineMVLHPerGene(Maps.newHashMap())
                .purple(TestPurpleInterpretationFactory.builder()
                        .germlineDrivers(Lists.newArrayList())
                        .allGermlineVariants(Lists.newArrayList())
                        .reportableGermlineVariants(Lists.newArrayList())
                        .additionalSuspectGermlineVariants(Lists.newArrayList())
                        .reportableGermlineFullLosses(Lists.newArrayList())
                        .build())
                .linx(TestLinxInterpretationFactory.builder()
                        .allGermlineStructuralVariants(Lists.newArrayList())
                        .allGermlineBreakends(Lists.newArrayList())
                        .reportableGermlineBreakends(Lists.newArrayList())
                        .germlineHomozygousDisruptions(Lists.newArrayList())
                        .build())
                .build();

        OrangeRecord converted = GermlineConversion.convertGermlineToSomatic(report);

        assertNull(converted.germlineMVLHPerGene());

        assertNull(converted.purple().germlineDrivers());
        assertNull(converted.purple().allGermlineVariants());
        assertNull(converted.purple().reportableGermlineVariants());
        assertNull(converted.purple().additionalSuspectGermlineVariants());
        assertNull(converted.purple().reportableGermlineFullLosses());

        assertNull(converted.linx().allGermlineStructuralVariants());
        assertNull(converted.linx().allGermlineBreakends());
        assertNull(converted.linx().reportableGermlineBreakends());
        assertNull(converted.linx().germlineHomozygousDisruptions());
    }

    @Test
    public void canConvertMinimalPurpleData()
    {
        OrangeRecord minimal = TestOrangeReportFactory.createMinimalTestReport();

        assertNotNull(GermlineConversion.convertPurpleGermline(true, minimal.purple()));
    }

    @Test
    public void canConvertPurple()
    {
        PurpleVariant somaticVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant reportableSomaticVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant suspectSomaticVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant germlineVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant reportableGermlineVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant suspectGermlineVariant = TestPurpleVariantFactory.builder().build();

        PurpleGainLoss somaticGainLoss = TestPurpleGainLossFactory.builder().build();
        PurpleGainLoss reportableSomaticGainLoss = TestPurpleGainLossFactory.builder().build();
        PurpleGainLoss reportableGermlineFullLoss = TestPurpleGainLossFactory.builder().build();
        PurpleLossOfHeterozygosity reportableGermlineLOH =
                TestPurpleLossOfHeterozygosityFactory.builder().gene(TEST_GENE1).minCopies(0.8).maxCopies(2.0).build();
        PurpleGeneCopyNumber geneCopyNumberForGermlineLOH =
                TestPurpleGeneCopyNumberFactory.builder().gene(TEST_GENE1).minCopyNumber(2.0).minMinorAlleleCopyNumber(0.9).maxCopyNumber(2.0).build();
        PurpleGeneCopyNumber suspectSomaticGeneCopyNumberWithLOH = TestPurpleGeneCopyNumberFactory.builder().gene(TEST_GENE2).build();

        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder().type(PurpleDriverType.AMP).build();
        PurpleDriver germlineMutationDriver = PurpleDriverTestFactory.builder().type(PurpleDriverType.GERMLINE_MUTATION).build();
        PurpleDriver germlineHomozygousDeletionDriver =
                PurpleDriverTestFactory.builder().type(PurpleDriverType.GERMLINE_DELETION).build();
        PurpleDriver germlineHeterozygousDeletionDriver =
                PurpleDriverTestFactory.builder().type(PurpleDriverType.GERMLINE_DELETION).gene(TEST_GENE1).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineMutationDriver, germlineHomozygousDeletionDriver, germlineHeterozygousDeletionDriver)
                .addAllSomaticVariants(somaticVariant, suspectSomaticVariant, reportableSomaticVariant)
                .addReportableSomaticVariants(reportableSomaticVariant)
                .addAdditionalSuspectSomaticVariants(suspectSomaticVariant)
                .addAllGermlineVariants(germlineVariant, suspectGermlineVariant, reportableGermlineVariant)
                .addReportableGermlineVariants(reportableGermlineVariant)
                .addAdditionalSuspectGermlineVariants(suspectGermlineVariant)
                .addAllSomaticGainsLosses(somaticGainLoss, reportableSomaticGainLoss)
                .addReportableSomaticGainsLosses(reportableSomaticGainLoss)
                .addReportableGermlineFullLosses(reportableGermlineFullLoss)
                .addSuspectGeneCopyNumbersWithLOH(suspectSomaticGeneCopyNumberWithLOH)
                .addReportableGermlineLossOfHeterozygosities(reportableGermlineLOH)
                .addAllSomaticGeneCopyNumbers(geneCopyNumberForGermlineLOH)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertTrue(converted.fit().qc().germlineAberrations().isEmpty());

        assertEquals(3, converted.somaticDrivers().size());
        assertNotNull(findByDriverType(converted.somaticDrivers(), PurpleDriverType.AMP));
        assertNotNull(findByDriverType(converted.somaticDrivers(), PurpleDriverType.MUTATION));
        assertNotNull(findByDriverType(converted.somaticDrivers(), PurpleDriverType.DEL));

        assertEquals(4, converted.allSomaticVariants().size());
        assertTrue(converted.allSomaticVariants().contains(reportableGermlineVariant));

        assertEquals(2, converted.reportableSomaticVariants().size());
        assertTrue(converted.reportableSomaticVariants().contains(reportableSomaticVariant));
        assertTrue(converted.reportableSomaticVariants().contains(reportableGermlineVariant));

        assertEquals(1, converted.additionalSuspectSomaticVariants().size());
        assertTrue(converted.additionalSuspectSomaticVariants().contains(suspectSomaticVariant));

        assertEquals(3, converted.allSomaticGainsLosses().size());
        assertTrue(converted.allSomaticGainsLosses().contains(reportableGermlineFullLoss));

        assertEquals(2, converted.reportableSomaticGainsLosses().size());
        assertTrue(converted.reportableSomaticGainsLosses().contains(reportableSomaticGainLoss));
        assertTrue(converted.reportableSomaticGainsLosses().contains(reportableGermlineFullLoss));

        assertEquals(2, converted.suspectGeneCopyNumbersWithLOH().size());
        assertTrue(converted.suspectGeneCopyNumbersWithLOH().contains(suspectSomaticGeneCopyNumberWithLOH));
        PurpleGeneCopyNumber convertedReportableGermlineLOH =
                TestPurpleGeneCopyNumberFactory.builder().gene(TEST_GENE1).minCopyNumber(0.8).minMinorAlleleCopyNumber(0.).maxCopyNumber(2.0).build();
        assertTrue(converted.suspectGeneCopyNumbersWithLOH().contains(convertedReportableGermlineLOH));

        PurpleRecord unreliableConverted = GermlineConversion.convertPurpleGermline(false, purple);
        assertEquals(1, unreliableConverted.somaticDrivers().size());
        assertNotNull(findByDriverType(unreliableConverted.somaticDrivers(), PurpleDriverType.AMP));

        assertEquals(1, unreliableConverted.reportableSomaticVariants().size());
        assertTrue(unreliableConverted.reportableSomaticVariants().contains(reportableSomaticVariant));

        assertEquals(1, unreliableConverted.reportableSomaticGainsLosses().size());
        assertTrue(unreliableConverted.reportableSomaticGainsLosses().contains(reportableSomaticGainLoss));
    }

    @Test
    public void mergesGermlineLOHIfAlreadySomaticLOH()
    {
        PurpleGeneCopyNumber suspectSomaticGeneCopyNumberWithLOH =
                TestPurpleGeneCopyNumberFactory.builder().gene(TEST_GENE1).minCopyNumber(0.7).minMinorAlleleCopyNumber(0.1).build();

        PurpleLossOfHeterozygosity germlineLOH = TestPurpleLossOfHeterozygosityFactory.builder().gene(TEST_GENE1).minCopies(0.8).build();
        PurpleGeneCopyNumber geneCopyNumber =
                TestPurpleGeneCopyNumberFactory.builder().gene(TEST_GENE1).minCopyNumber(2.0).minMinorAlleleCopyNumber(0.9).build();
        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder().type(PurpleDriverType.GERMLINE_DELETION).gene(TEST_GENE1).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses()
                .addSuspectGeneCopyNumbersWithLOH(suspectSomaticGeneCopyNumberWithLOH)
                .addReportableGermlineLossOfHeterozygosities(germlineLOH)
                .addAllSomaticGeneCopyNumbers(geneCopyNumber)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(0, converted.somaticDrivers().size());

        assertEquals(1, converted.suspectGeneCopyNumbersWithLOH().size());
        assertEquals(0.7, converted.suspectGeneCopyNumbersWithLOH().get(0).minCopyNumber(), EPSILON);
        assertEquals(0., converted.suspectGeneCopyNumbersWithLOH().get(0).minMinorAlleleCopyNumber(), EPSILON);
    }

    @Test
    public void doesNotConvertGermlineLOHIfAlsoGermlineFullLoss()
    {
        PurpleGainLoss reportableGermlineFullLoss = TestPurpleGainLossFactory.builder().gene(TEST_GENE1).build();

        PurpleLossOfHeterozygosity germlineLOH = TestPurpleLossOfHeterozygosityFactory.builder().gene(TEST_GENE1).minCopies(0.8).build();
        PurpleGeneCopyNumber geneCopyNumber =
                TestPurpleGeneCopyNumberFactory.builder().gene(TEST_GENE1).minCopyNumber(2.0).minMinorAlleleCopyNumber(0.9).build();
        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder().type(PurpleDriverType.GERMLINE_DELETION).gene(TEST_GENE1).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses()
                .addReportableGermlineFullLosses(reportableGermlineFullLoss)
                .addReportableGermlineLossOfHeterozygosities(germlineLOH)
                .addAllSomaticGeneCopyNumbers(geneCopyNumber)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(1, converted.somaticDrivers().size());
        PurpleDriver convertedDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .gene(TEST_GENE1)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .build();
        assertTrue(converted.somaticDrivers().contains(convertedDriver));

        assertEquals(0, converted.suspectGeneCopyNumbersWithLOH().size());
    }

    @Test
    public void doesNotConvertGermlineLOHIfAlsoFullLossInTumor()
    {
        PurpleLossOfHeterozygosity germlineLOH = TestPurpleLossOfHeterozygosityFactory.builder().gene(TEST_GENE1).minCopies(0.8).build();
        PurpleGeneCopyNumber geneCopyNumber =
                TestPurpleGeneCopyNumberFactory.builder().gene(TEST_GENE1).minCopyNumber(0.3).minMinorAlleleCopyNumber(0.).build();
        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder().type(PurpleDriverType.GERMLINE_DELETION).gene(TEST_GENE1).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses()
                .addReportableGermlineLossOfHeterozygosities(germlineLOH)
                .addAllSomaticGeneCopyNumbers(geneCopyNumber)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(0, converted.somaticDrivers().size());
        assertEquals(0, converted.suspectGeneCopyNumbersWithLOH().size());
    }

    @Test
    public void convertsMultipleHeterozygousGermlineDeletionsCorrectly()
    {
        PurpleLossOfHeterozygosity germlineLOH1 = TestPurpleLossOfHeterozygosityFactory.builder().gene(TEST_GENE1).minCopies(0.8).build();
        PurpleLossOfHeterozygosity germlineLOH2 = TestPurpleLossOfHeterozygosityFactory.builder().gene(TEST_GENE1).minCopies(0.9).build();
        PurpleGeneCopyNumber geneCopyNumber =
                TestPurpleGeneCopyNumberFactory.builder().gene(TEST_GENE1).minCopyNumber(2.0).minMinorAlleleCopyNumber(0.9).build();
        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder().type(PurpleDriverType.GERMLINE_DELETION).gene(TEST_GENE1).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses()
                .addReportableGermlineLossOfHeterozygosities(germlineLOH1, germlineLOH2)
                .addAllSomaticGeneCopyNumbers(geneCopyNumber)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(0, converted.somaticDrivers().size());

        assertEquals(1, converted.suspectGeneCopyNumbersWithLOH().size());
        PurpleGeneCopyNumber convertedGeneCopyNumber =
                TestPurpleGeneCopyNumberFactory.builder().gene(TEST_GENE1).minCopyNumber(0.8).minMinorAlleleCopyNumber(0.).build();
        assertTrue(converted.suspectGeneCopyNumbersWithLOH().contains(convertedGeneCopyNumber));
    }

    @Test
    public void mergesSomaticAndGermlinePartialLosses()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss somaticLoss = TestPurpleGainLossFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.PARTIAL_LOSS)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.1)
                .maxCopies(1.0)
                .build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss germlineLoss = TestPurpleGainLossFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.PARTIAL_LOSS)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.2)
                .maxCopies(0.9)
                .build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses(somaticLoss)
                .addReportableGermlineFullLosses(germlineLoss)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(1, converted.somaticDrivers().size());
        assertEquals(somaticDriver, converted.somaticDrivers().get(0));

        assertEquals(1, converted.reportableSomaticGainsLosses().size());
        assertEquals(0.1, converted.reportableSomaticGainsLosses().get(0).minCopies(), EPSILON);
        assertEquals(1.0, converted.reportableSomaticGainsLosses().get(0).maxCopies(), EPSILON);
        assertEquals(CopyNumberInterpretation.PARTIAL_LOSS, converted.reportableSomaticGainsLosses().get(0).interpretation());
    }

    @Test
    public void mergesSomaticAndGermlineFullLosses()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss somaticLoss = TestPurpleGainLossFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.1)
                .maxCopies(0.3)
                .build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss germlineLoss = TestPurpleGainLossFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.2)
                .maxCopies(0.4)
                .build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses(somaticLoss)
                .addReportableGermlineFullLosses(germlineLoss)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(1, converted.somaticDrivers().size());
        assertEquals(somaticDriver, converted.somaticDrivers().get(0));

        assertEquals(1, converted.reportableSomaticGainsLosses().size());
        assertEquals(0.1, converted.reportableSomaticGainsLosses().get(0).minCopies(), EPSILON);
        assertEquals(0.4, converted.reportableSomaticGainsLosses().get(0).maxCopies(), EPSILON);
        assertEquals(CopyNumberInterpretation.FULL_LOSS, converted.reportableSomaticGainsLosses().get(0).interpretation());
    }

    @Test
    public void mergesSomaticPartialAndGermlineFullLoss()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss somaticLoss = TestPurpleGainLossFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.PARTIAL_LOSS)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.1)
                .maxCopies(1.0)
                .build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss germlineLoss = TestPurpleGainLossFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.2)
                .maxCopies(0.4)
                .build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses(somaticLoss)
                .addReportableGermlineFullLosses(germlineLoss)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(1, converted.somaticDrivers().size());
        assertEquals(somaticDriver, converted.somaticDrivers().get(0));

        assertEquals(1, converted.reportableSomaticGainsLosses().size());
        assertEquals(0.1, converted.reportableSomaticGainsLosses().get(0).minCopies(), EPSILON);
        assertEquals(0.4, converted.reportableSomaticGainsLosses().get(0).maxCopies(), EPSILON);
        assertEquals(CopyNumberInterpretation.FULL_LOSS, converted.reportableSomaticGainsLosses().get(0).interpretation());
    }

    @Test
    public void mergesSomaticFullAndGermlinePartialLoss()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss somaticLoss = TestPurpleGainLossFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.FULL_LOSS)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.1)
                .maxCopies(0.3)
                .build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss germlineLoss = TestPurpleGainLossFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.PARTIAL_LOSS)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.2)
                .maxCopies(1.1)
                .build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses(somaticLoss)
                .addReportableGermlineFullLosses(germlineLoss)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(1, converted.somaticDrivers().size());
        assertEquals(somaticDriver, converted.somaticDrivers().get(0));

        assertEquals(1, converted.reportableSomaticGainsLosses().size());
        assertEquals(0.1, converted.reportableSomaticGainsLosses().get(0).minCopies(), EPSILON);
        assertEquals(0.3, converted.reportableSomaticGainsLosses().get(0).maxCopies(), EPSILON);
        assertEquals(CopyNumberInterpretation.FULL_LOSS, converted.reportableSomaticGainsLosses().get(0).interpretation());
    }

    @Test
    public void doesNotMergeSomaticAndGermlineLossesOnDifferentTranscript()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .gene(TEST_GENE1)
                .transcript(TEST_TRANSCRIPT1)
                .isCanonical(false)
                .build();
        PurpleGainLoss somaticLoss =
                TestPurpleGainLossFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT1).minCopies(0.1).maxCopies(1.0).build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .transcript(TEST_TRANSCRIPT2)
                .isCanonical(true)
                .build();

        PurpleGainLoss germlineLoss =
                TestPurpleGainLossFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT2).minCopies(0.2).maxCopies(0.9).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses(somaticLoss)
                .addReportableGermlineFullLosses(germlineLoss)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(2, converted.somaticDrivers().size());
        assertEquals(somaticDriver, converted.somaticDrivers().get(0));

        assertEquals(2, converted.reportableSomaticGainsLosses().size());
    }

    @Test
    public void doesNotMergeSomaticAmpWithGermlineLoss()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.AMP)
                .likelihoodMethod(PurpleLikelihoodMethod.AMP)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss somaticAmp =
                TestPurpleGainLossFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT1).minCopies(20.).maxCopies(25.).build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss germlineLoss =
                TestPurpleGainLossFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT1).minCopies(0.2).maxCopies(0.9).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses(somaticAmp)
                .addReportableGermlineFullLosses(germlineLoss)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(2, converted.somaticDrivers().size());
        assertEquals(2, converted.reportableSomaticGainsLosses().size());
    }

    @Test
    public void doesNotMergeSomaticPartialAmpWithGermlineLoss()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.PARTIAL_AMP)
                .likelihoodMethod(PurpleLikelihoodMethod.AMP)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss somaticPartialAmp =
                TestPurpleGainLossFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT1).minCopies(2.0).maxCopies(25.).build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainLoss germlineLoss =
                TestPurpleGainLossFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT1).minCopies(0.2).maxCopies(0.3).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addReportableSomaticGainsLosses(somaticPartialAmp)
                .addReportableGermlineFullLosses(germlineLoss)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(2, converted.somaticDrivers().size());
        assertEquals(2, converted.reportableSomaticGainsLosses().size());
    }

    @NotNull
    private static PurpleFit createWithGermlineAberration()
    {
        PurpleFit base = TestPurpleInterpretationFactory.createMinimalTestFitData();
        return ImmutablePurpleFit.builder()
                .from(base)
                .qc(TestPurpleQCFactory.builder().from(base.qc()).addGermlineAberrations(PurpleGermlineAberration.MOSAIC_X).build())
                .build();
    }

    @Nullable
    private static PurpleDriver findByDriverType(@NotNull List<PurpleDriver> drivers, @NotNull PurpleDriverType driverTypeToFind)
    {
        for(PurpleDriver driver : drivers)
        {
            if(driver.type() == driverTypeToFind)
            {
                return driver;
            }
        }

        return null;
    }

    @Test
    public void canMergeMutationDrivers()
    {
        PurpleDriver somaticDriver1 = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.MUTATION)
                .gene("gene 1")
                .transcript("transcript 1")
                .driverLikelihood(0.5)
                .likelihoodMethod(PurpleLikelihoodMethod.DNDS)
                .build();

        PurpleDriver somaticDriver2 =
                PurpleDriverTestFactory.builder().type(PurpleDriverType.MUTATION).gene("gene 1").transcript("transcript 2").build();

        PurpleDriver germlineDriver1 = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_MUTATION)
                .gene("gene 1")
                .transcript("transcript 1")
                .driverLikelihood(0.8)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .build();

        PurpleDriver germlineDriver2 = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_MUTATION)
                .gene("gene 2")
                .transcript("transcript 1")
                .build();

        List<PurpleDriver> merged = GermlineConversion.mergeGermlineDriversIntoSomatic(
                Lists.newArrayList(somaticDriver1, somaticDriver2),
                Lists.newArrayList(germlineDriver1, germlineDriver2),
                Lists.newArrayList()
        );

        assertEquals(3, merged.size());
        PurpleDriver driver1 = findByGeneTranscript(merged, "gene 1", "transcript 1");
        assertNotNull(driver1);
        assertEquals(PurpleDriverType.MUTATION, driver1.type());
        assertEquals(PurpleLikelihoodMethod.HOTSPOT, driver1.likelihoodMethod());
        assertEquals(0.8, driver1.driverLikelihood(), EPSILON);

        PurpleDriver driver2 = findByGeneTranscript(merged, "gene 1", "transcript 2");
        assertEquals(somaticDriver2, driver2);

        PurpleDriver driver3 = findByGeneTranscript(merged, "gene 2", "transcript 1");
        assertEquals(PurpleDriverType.MUTATION, driver3.type());
        assertEquals(PurpleLikelihoodMethod.HOTSPOT, driver3.likelihoodMethod());
    }

    @Test
    public void canConvertAllTypesOfGermlineDrivers()
    {
        PurpleDriver somaticDriver =
                PurpleDriverTestFactory.builder().gene("gene 1").transcript("transcript 1").type(PurpleDriverType.DEL).build();

        List<PurpleDriver> mergedNoGermline =
                GermlineConversion.mergeGermlineDriversIntoSomatic(Lists.newArrayList(somaticDriver), null, Lists.newArrayList());
        assertEquals(1, mergedNoGermline.size());
        assertTrue(mergedNoGermline.contains(somaticDriver));

        PurpleDriver germlineDriver1 = PurpleDriverTestFactory.builder()
                .gene("gene 1")
                .transcript("transcript 1")
                .type(PurpleDriverType.GERMLINE_DELETION)
                .build();

        PurpleDriver germlineDriver2 = PurpleDriverTestFactory.builder()
                .gene("gene 2")
                .transcript("transcript 2")
                .type(PurpleDriverType.GERMLINE_DELETION)
                .build();
        PurpleGainLoss fullGermlineLossForDriver2 =
                TestPurpleGainLossFactory.builder().gene("gene 2").interpretation(CopyNumberInterpretation.FULL_LOSS).build();

        PurpleDriver germlineDriver3 = PurpleDriverTestFactory.builder()
                .gene("gene 3")
                .transcript("transcript 3")
                .type(PurpleDriverType.GERMLINE_MUTATION)
                .build();

        List<PurpleDriver> merged = GermlineConversion.mergeGermlineDriversIntoSomatic(
                Lists.newArrayList(somaticDriver),
                Lists.newArrayList(germlineDriver1, germlineDriver2, germlineDriver3),
                Lists.newArrayList(fullGermlineLossForDriver2)
        );

        assertEquals(3, merged.size());
        assertTrue(mergedNoGermline.contains(somaticDriver));

        PurpleDriver germlineDeletionDriver = findByGeneTranscript(merged, germlineDriver2.gene(), germlineDriver2.transcript());
        assertEquals(PurpleDriverType.DEL, germlineDeletionDriver.type());

        PurpleDriver germlineMutationDriver = findByGeneTranscript(merged, germlineDriver3.gene(), germlineDriver3.transcript());
        assertEquals(PurpleDriverType.MUTATION, germlineMutationDriver.type());
    }

    @Nullable
    private static PurpleDriver findByGeneTranscript(@NotNull List<PurpleDriver> drivers, @NotNull String geneToFind,
            @NotNull String transcriptToFind)
    {
        for(PurpleDriver driver : drivers)
        {
            if(driver.gene().equals(geneToFind) && driver.transcript().equals(transcriptToFind))
            {
                return driver;
            }
        }
        return null;
    }

    @Test
    public void canConvertLinx()
    {
        LinxSvAnnotation somaticStructuralVariant1 = LinxOrangeTestFactory.svAnnotationBuilder().svId(1).clusterId(5).build();
        LinxSvAnnotation somaticStructuralVariant2 = LinxOrangeTestFactory.svAnnotationBuilder().svId(2).clusterId(6).build();
        LinxBreakend somaticBreakend = LinxOrangeTestFactory.breakendBuilder().id(8).svId(1).build();
        LinxBreakend reportableSomaticBreakend = LinxOrangeTestFactory.breakendBuilder().id(9).svId(2).build();

        LinxSvAnnotation germlineStructuralVariant1 = LinxOrangeTestFactory.svAnnotationBuilder().svId(1).clusterId(5).build();
        LinxSvAnnotation germlineStructuralVariant2 = LinxOrangeTestFactory.svAnnotationBuilder().svId(2).clusterId(6).build();
        LinxSvAnnotation germlineStructuralVariant3 = LinxOrangeTestFactory.svAnnotationBuilder().svId(3).clusterId(6).build();
        LinxBreakend germlineBreakend = LinxOrangeTestFactory.breakendBuilder().id(8).svId(1).build();
        LinxBreakend reportableGermlineBreakend = LinxOrangeTestFactory.breakendBuilder().id(9).svId(2).build();

        LinxHomozygousDisruption germlineHomozygousDisruption = LinxOrangeTestFactory.homozygousDisruptionBuilder().build();

        LinxRecord linx = TestLinxInterpretationFactory.builder()
                .addAllSomaticStructuralVariants(somaticStructuralVariant1, somaticStructuralVariant2)
                .addAllSomaticBreakends(somaticBreakend, reportableSomaticBreakend)
                .addReportableSomaticBreakends(reportableSomaticBreakend)
                .addAllGermlineStructuralVariants(germlineStructuralVariant1, germlineStructuralVariant2, germlineStructuralVariant3)
                .addAllGermlineBreakends(germlineBreakend, reportableGermlineBreakend)
                .addReportableGermlineBreakends(reportableGermlineBreakend)
                .addGermlineHomozygousDisruptions(germlineHomozygousDisruption)
                .build();

        LinxRecord converted = GermlineConversion.convertLinxGermline(true, linx);
        assertEquals(5, converted.allSomaticStructuralVariants().size());
        assertEquals(5, GermlineConversion.findMaxSvId(converted.allSomaticStructuralVariants()));
        assertEquals(8, GermlineConversion.findMaxClusterId(converted.allSomaticStructuralVariants()));

        assertEquals(3, converted.allSomaticBreakends().size());
        assertEquals(11, GermlineConversion.findMaxBreakendId(converted.allSomaticBreakends()));

        assertEquals(2, converted.reportableSomaticBreakends().size());
        assertEquals(11, GermlineConversion.findMaxBreakendId(converted.reportableSomaticBreakends()));

        assertEquals(1, converted.somaticHomozygousDisruptions().size());
        assertTrue(converted.somaticHomozygousDisruptions().contains(germlineHomozygousDisruption));
    }

    @Test
    public void shouldAdjustClonalLikelihoodWhenConvertingVariantsToSomatic()
    {
        PurpleVariant germlineVariant1 = TestPurpleVariantFactory.builder().variantCopyNumber(0).build();
        final List<PurpleVariant> somaticVariants1 = GermlineConversion.toSomaticVariants(Lists.newArrayList(germlineVariant1));
        assertEquals(1, somaticVariants1.size());
        assertEquals(1.0, somaticVariants1.get(0).subclonalLikelihood(), EPSILON);

        PurpleVariant germlineVariant2 = TestPurpleVariantFactory.builder().variantCopyNumber(1).build();
        final List<PurpleVariant> somaticVariants2 = GermlineConversion.toSomaticVariants(Lists.newArrayList(germlineVariant2));
        assertEquals(1, somaticVariants1.size());
        assertEquals(0.0, somaticVariants2.get(0).subclonalLikelihood(), EPSILON);
    }
}