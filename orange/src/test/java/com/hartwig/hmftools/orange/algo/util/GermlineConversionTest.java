package com.hartwig.hmftools.orange.algo.util;

import static com.hartwig.hmftools.orange.algo.linx.HomozygousDisruptionFactoryTest.linxHomozygousDisruptionBuilder;
import static com.hartwig.hmftools.orange.algo.purple.TumorStatsFactoryTest.createMinimalTumorStatsBuilder;

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
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGermlineAberration;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.TumorStats;
import com.hartwig.hmftools.orange.TestOrangeReportFactory;
import com.hartwig.hmftools.orange.algo.linx.LinxOrangeTestFactory;
import com.hartwig.hmftools.orange.algo.linx.TestLinxInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGainDeletionFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGeneCopyNumberFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleLossOfHeterozygosityFactory;
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
                        .otherGermlineVariants(Lists.newArrayList())
                        .driverGermlineVariants(Lists.newArrayList())
                        .build())
                .linx(TestLinxInterpretationFactory.builder()
                        .allGermlineStructuralVariants(Lists.newArrayList())
                        .otherGermlineBreakends(Lists.newArrayList())
                        .driverGermlineBreakends(Lists.newArrayList())
                        .germlineHomozygousDisruptions(Lists.newArrayList())
                        .build())
                .build();

        OrangeRecord converted = GermlineConversion.convertGermlineToSomatic(report);

        assertNull(converted.germlineMVLHPerGene());

        assertNull(converted.purple().germlineDrivers());
        assertNull(converted.purple().otherGermlineVariants());
        assertNull(converted.purple().driverGermlineVariants());
        assertNull(converted.purple().driverGermlineDeletions());

        assertNull(converted.linx().allGermlineStructuralVariants());
        assertNull(converted.linx().otherGermlineBreakends());
        assertNull(converted.linx().driverGermlineBreakends());
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
        PurpleVariant reportableSomaticVariant = TestPurpleVariantFactory.builder()
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().reported(true).build())
                .build();
        PurpleVariant suspectSomaticVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant germlineVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant reportableGermlineVariant = TestPurpleVariantFactory.builder()
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().reported(true).build())
                .build();
        PurpleVariant suspectGermlineVariant = TestPurpleVariantFactory.builder().build();

        PurpleGainDeletion reportableSomaticGainDel = TestPurpleGainDeletionFactory.builder().build();
        PurpleGainDeletion reportableGermlineFullDel = TestPurpleGainDeletionFactory.builder().build();
        PurpleLossOfHeterozygosity reportableGermlineLOH =
                TestPurpleLossOfHeterozygosityFactory.builder().gene(TEST_GENE1).minCopies(0.8).maxCopies(2.0).build();
        PurpleGeneCopyNumber geneCopyNumberForGermlineLOH =
                TestPurpleGeneCopyNumberFactory.builder()
                        .gene(TEST_GENE1)
                        .minCopyNumber(2.0)
                        .minMinorAlleleCopyNumber(0.9)
                        .maxCopyNumber(2.0)
                        .build();

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
                .addOtherSomaticVariants(somaticVariant, suspectSomaticVariant, reportableSomaticVariant)
                .addDriverSomaticVariants(reportableSomaticVariant)
                .addOtherGermlineVariants(germlineVariant, suspectGermlineVariant, reportableGermlineVariant)
                .addDriverGermlineVariants(reportableGermlineVariant)
                .addDriverSomaticGainsDels(reportableSomaticGainDel)
                .addDriverGermlineDeletions(reportableGermlineFullDel)
                .addDriverGermlineLossOfHeterozygosities(reportableGermlineLOH)
                .addSomaticGeneCopyNumbers(geneCopyNumberForGermlineLOH)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertTrue(converted.fit().qc().germlineAberrations().isEmpty());

        assertEquals(3, converted.somaticDrivers().size());
        assertNotNull(findByDriverType(converted.somaticDrivers(), PurpleDriverType.AMP));
        assertNotNull(findByDriverType(converted.somaticDrivers(), PurpleDriverType.MUTATION));
        assertNotNull(findByDriverType(converted.somaticDrivers(), PurpleDriverType.DEL));

        assertEquals(4, converted.otherSomaticVariants().size());
        assertTrue(converted.otherSomaticVariants().contains(reportableGermlineVariant));

        assertEquals(2, converted.driverSomaticVariants().size());
        assertTrue(converted.driverSomaticVariants().contains(reportableSomaticVariant));
        assertTrue(converted.driverSomaticVariants().contains(reportableGermlineVariant));

        assertEquals(2, converted.driverSomaticGainsDels().size());
        assertTrue(converted.driverSomaticGainsDels().contains(reportableSomaticGainDel));
        assertTrue(converted.driverSomaticGainsDels().contains(reportableGermlineFullDel));

        PurpleGeneCopyNumber convertedReportableGermlineLOH =
                TestPurpleGeneCopyNumberFactory.builder()
                        .gene(TEST_GENE1)
                        .minCopyNumber(0.8)
                        .minMinorAlleleCopyNumber(0.)
                        .maxCopyNumber(2.0)
                        .build();
        PurpleRecord unreliableConverted = GermlineConversion.convertPurpleGermline(false, purple);
        assertEquals(1, unreliableConverted.somaticDrivers().size());
        assertNotNull(findByDriverType(unreliableConverted.somaticDrivers(), PurpleDriverType.AMP));

        assertEquals(1, unreliableConverted.driverSomaticVariants().size());
        assertTrue(unreliableConverted.driverSomaticVariants().contains(reportableSomaticVariant));

        assertEquals(1, unreliableConverted.driverSomaticGainsDels().size());
        assertTrue(unreliableConverted.driverSomaticGainsDels().contains(reportableSomaticGainDel));
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
                .addDriverSomaticGainsDels()
                .addDriverGermlineLossOfHeterozygosities(germlineLOH)
                .addSomaticGeneCopyNumbers(geneCopyNumber)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(0, converted.somaticDrivers().size());
    }

    @Test
    public void doesNotConvertGermlineLOHIfAlsoGermlineFullLoss()
    {
        PurpleGainDeletion reportableGermlineFullDel = TestPurpleGainDeletionFactory.builder().gene(TEST_GENE1).build();

        PurpleLossOfHeterozygosity germlineLOH = TestPurpleLossOfHeterozygosityFactory.builder().gene(TEST_GENE1).minCopies(0.8).build();
        PurpleGeneCopyNumber geneCopyNumber =
                TestPurpleGeneCopyNumberFactory.builder().gene(TEST_GENE1).minCopyNumber(2.0).minMinorAlleleCopyNumber(0.9).build();
        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder().type(PurpleDriverType.GERMLINE_DELETION).gene(TEST_GENE1).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addGermlineDrivers(germlineDriver)
                .addDriverSomaticGainsDels()
                .addDriverGermlineDeletions(reportableGermlineFullDel)
                .addDriverGermlineLossOfHeterozygosities(germlineLOH)
                .addSomaticGeneCopyNumbers(geneCopyNumber)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(1, converted.somaticDrivers().size());
        PurpleDriver convertedDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .gene(TEST_GENE1)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .build();
        assertTrue(converted.somaticDrivers().contains(convertedDriver));
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
                .addDriverSomaticGainsDels()
                .addDriverGermlineLossOfHeterozygosities(germlineLOH)
                .addSomaticGeneCopyNumbers(geneCopyNumber)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(0, converted.somaticDrivers().size());
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
                .addDriverSomaticGainsDels()
                .addDriverGermlineLossOfHeterozygosities(germlineLOH1, germlineLOH2)
                .addSomaticGeneCopyNumbers(geneCopyNumber)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(0, converted.somaticDrivers().size());
    }

    @Test
    public void mergesSomaticAndGermlinePartialDeletions()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion somaticDel = TestPurpleGainDeletionFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.PARTIAL_DEL)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.1)
                .maxCopies(1.0)
                .build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion germlineDel = TestPurpleGainDeletionFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.PARTIAL_DEL)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.2)
                .maxCopies(0.9)
                .build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addDriverSomaticGainsDels(somaticDel)
                .addDriverGermlineDeletions(germlineDel)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(1, converted.somaticDrivers().size());
        assertEquals(somaticDriver, converted.somaticDrivers().get(0));

        assertEquals(1, converted.driverSomaticGainsDels().size());
        assertEquals(0.1, converted.driverSomaticGainsDels().get(0).minCopies(), EPSILON);
        assertEquals(1.0, converted.driverSomaticGainsDels().get(0).maxCopies(), EPSILON);
        assertEquals(CopyNumberInterpretation.PARTIAL_DEL, converted.driverSomaticGainsDels().get(0).interpretation());
    }

    @Test
    public void mergesSomaticAndGermlineFullDeletions()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion somaticDel = TestPurpleGainDeletionFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.FULL_DEL)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.1)
                .maxCopies(0.3)
                .build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion germlineDel = TestPurpleGainDeletionFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.FULL_DEL)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.2)
                .maxCopies(0.4)
                .build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addDriverSomaticGainsDels(somaticDel)
                .addDriverGermlineDeletions(germlineDel)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(1, converted.somaticDrivers().size());
        assertEquals(somaticDriver, converted.somaticDrivers().get(0));

        assertEquals(1, converted.driverSomaticGainsDels().size());
        assertEquals(0.1, converted.driverSomaticGainsDels().get(0).minCopies(), EPSILON);
        assertEquals(0.4, converted.driverSomaticGainsDels().get(0).maxCopies(), EPSILON);
        assertEquals(CopyNumberInterpretation.FULL_DEL, converted.driverSomaticGainsDels().get(0).interpretation());
    }

    @Test
    public void mergesSomaticPartialAndGermlineFullDeletion()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion somaticDel = TestPurpleGainDeletionFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.PARTIAL_DEL)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.1)
                .maxCopies(1.0)
                .build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion germlineDel = TestPurpleGainDeletionFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.FULL_DEL)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.2)
                .maxCopies(0.4)
                .build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addDriverSomaticGainsDels(somaticDel)
                .addDriverGermlineDeletions(germlineDel)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(1, converted.somaticDrivers().size());
        assertEquals(somaticDriver, converted.somaticDrivers().get(0));

        assertEquals(1, converted.driverSomaticGainsDels().size());
        assertEquals(0.1, converted.driverSomaticGainsDels().get(0).minCopies(), EPSILON);
        assertEquals(0.4, converted.driverSomaticGainsDels().get(0).maxCopies(), EPSILON);
        assertEquals(CopyNumberInterpretation.FULL_DEL, converted.driverSomaticGainsDels().get(0).interpretation());
    }

    @Test
    public void mergesSomaticFullAndGermlinePartialDeletion()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion somaticDel = TestPurpleGainDeletionFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.FULL_DEL)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.1)
                .maxCopies(0.3)
                .build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion germlineDel = TestPurpleGainDeletionFactory.builder()
                .gene(TEST_GENE1)
                .interpretation(CopyNumberInterpretation.PARTIAL_DEL)
                .transcript(TEST_TRANSCRIPT1)
                .minCopies(0.2)
                .maxCopies(1.1)
                .build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addDriverSomaticGainsDels(somaticDel)
                .addDriverGermlineDeletions(germlineDel)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(1, converted.somaticDrivers().size());
        assertEquals(somaticDriver, converted.somaticDrivers().get(0));

        assertEquals(1, converted.driverSomaticGainsDels().size());
        assertEquals(0.1, converted.driverSomaticGainsDels().get(0).minCopies(), EPSILON);
        assertEquals(0.3, converted.driverSomaticGainsDels().get(0).maxCopies(), EPSILON);
        assertEquals(CopyNumberInterpretation.FULL_DEL, converted.driverSomaticGainsDels().get(0).interpretation());
    }

    @Test
    public void doesNotMergeSomaticAndGermlineDeletionsOnDifferentTranscript()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.DEL)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .gene(TEST_GENE1)
                .transcript(TEST_TRANSCRIPT1)
                .isCanonical(false)
                .build();
        PurpleGainDeletion somaticDel =
                TestPurpleGainDeletionFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT1).minCopies(0.1).maxCopies(1.0).build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .transcript(TEST_TRANSCRIPT2)
                .isCanonical(true)
                .build();

        PurpleGainDeletion germlineDel =
                TestPurpleGainDeletionFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT2).minCopies(0.2).maxCopies(0.9).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addDriverSomaticGainsDels(somaticDel)
                .addDriverGermlineDeletions(germlineDel)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(2, converted.somaticDrivers().size());
        assertEquals(somaticDriver, converted.somaticDrivers().get(0));

        assertEquals(2, converted.driverSomaticGainsDels().size());
    }

    @Test
    public void doesNotMergeSomaticAmpWithGermlineDeletion()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.AMP)
                .likelihoodMethod(PurpleLikelihoodMethod.AMP)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion somaticAmp =
                TestPurpleGainDeletionFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT1).minCopies(20.).maxCopies(25.).build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion germlineDel =
                TestPurpleGainDeletionFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT1).minCopies(0.2).maxCopies(0.9).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addDriverSomaticGainsDels(somaticAmp)
                .addDriverGermlineDeletions(germlineDel)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(2, converted.somaticDrivers().size());
        assertEquals(2, converted.driverSomaticGainsDels().size());
    }

    @Test
    public void doesNotMergeSomaticPartialAmpWithGermlineDeletion()
    {
        PurpleDriver somaticDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.PARTIAL_AMP)
                .likelihoodMethod(PurpleLikelihoodMethod.AMP)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion somaticPartialAmp =
                TestPurpleGainDeletionFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT1).minCopies(2.0).maxCopies(25.).build();

        PurpleDriver germlineDriver = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.GERMLINE_DELETION)
                .likelihoodMethod(PurpleLikelihoodMethod.GERMLINE)
                .gene(TEST_GENE1)
                .build();
        PurpleGainDeletion germlineDel =
                TestPurpleGainDeletionFactory.builder().gene(TEST_GENE1).transcript(TEST_TRANSCRIPT1).minCopies(0.2).maxCopies(0.3).build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineDriver)
                .addDriverSomaticGainsDels(somaticPartialAmp)
                .addDriverGermlineDeletions(germlineDel)
                .build();

        PurpleRecord converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertEquals(2, converted.somaticDrivers().size());
        assertEquals(2, converted.driverSomaticGainsDels().size());
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
        PurpleGainDeletion fullgermlineDelForDriver2 =
                TestPurpleGainDeletionFactory.builder().gene("gene 2").interpretation(CopyNumberInterpretation.FULL_DEL).build();

        PurpleDriver germlineDriver3 = PurpleDriverTestFactory.builder()
                .gene("gene 3")
                .transcript("transcript 3")
                .type(PurpleDriverType.GERMLINE_MUTATION)
                .build();

        List<PurpleDriver> merged = GermlineConversion.mergeGermlineDriversIntoSomatic(
                Lists.newArrayList(somaticDriver),
                Lists.newArrayList(germlineDriver1, germlineDriver2, germlineDriver3),
                Lists.newArrayList(fullgermlineDelForDriver2)
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
        LinxBreakend reportableSomaticBreakend = LinxOrangeTestFactory.breakendBuilder().id(9).svId(2).reported(true).build();

        LinxSvAnnotation germlineStructuralVariant1 = LinxOrangeTestFactory.svAnnotationBuilder().svId(1).clusterId(5).build();
        LinxSvAnnotation germlineStructuralVariant2 = LinxOrangeTestFactory.svAnnotationBuilder().svId(2).clusterId(6).build();
        LinxSvAnnotation germlineStructuralVariant3 = LinxOrangeTestFactory.svAnnotationBuilder().svId(3).clusterId(6).build();
        LinxBreakend germlineBreakend = LinxOrangeTestFactory.breakendBuilder().id(8).svId(1).build();
        LinxBreakend reportableGermlineBreakend = LinxOrangeTestFactory.breakendBuilder().id(9).svId(2).reported(true).build();

        LinxHomozygousDisruption germlineHomozygousDisruption = linxHomozygousDisruptionBuilder().build();

        LinxRecord linx = TestLinxInterpretationFactory.builder()
                .addAllSomaticStructuralVariants(somaticStructuralVariant1, somaticStructuralVariant2)
                .addOtherSomaticBreakends(somaticBreakend, reportableSomaticBreakend)
                .addDriverSomaticBreakends(reportableSomaticBreakend)
                .addAllGermlineStructuralVariants(germlineStructuralVariant1, germlineStructuralVariant2, germlineStructuralVariant3)
                .addOtherGermlineBreakends(germlineBreakend, reportableGermlineBreakend)
                .addDriverGermlineBreakends(reportableGermlineBreakend)
                .addGermlineHomozygousDisruptions(germlineHomozygousDisruption)
                .build();

        LinxRecord converted = GermlineConversion.convertLinxGermline(true, linx);
        assertEquals(5, converted.allSomaticStructuralVariants().size());
        assertEquals(5, GermlineConversion.findMaxSvId(converted.allSomaticStructuralVariants()));
        assertEquals(8, GermlineConversion.findMaxClusterId(converted.allSomaticStructuralVariants()));

        assertEquals(3, converted.otherSomaticBreakends().size());
        assertEquals(11, GermlineConversion.findMaxBreakendId(converted.otherSomaticBreakends()));

        assertEquals(2, converted.driverSomaticBreakends().size());
        assertEquals(11, GermlineConversion.findMaxBreakendId(converted.driverSomaticBreakends()));

        assertEquals(1, converted.somaticHomozygousDisruptions().size());
        assertTrue(converted.somaticHomozygousDisruptions().contains(germlineHomozygousDisruption));
    }

    @Test
    public void canAdjustClonalLikelihoodWhenConvertingVariantsToSomatic()
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

    @Test
    public void canMergeTumorStats()
    {
        PurpleVariant somaticVariant = TestPurpleVariantFactory.builder()
                .hotspot(HotspotType.HOTSPOT)
                .build();

        PurpleVariant reportableGermlineVariant = TestPurpleVariantFactory.builder()
                .hotspot(HotspotType.HOTSPOT)
                .canonicalImpact(TestPurpleVariantFactory.impactBuilder().reported(true).build())
                .build();

        PurpleRecord purple = TestPurpleInterpretationFactory.builder()
                .addOtherSomaticVariants(somaticVariant)
                .addDriverGermlineVariants(reportableGermlineVariant)
                .tumorStats(createMinimalTumorStatsBuilder().hotspotMutationCount(1).build())
                .build();

        TumorStats mergedTumorStats = GermlineConversion.mergeTumorStats(purple);
        assertEquals(2, mergedTumorStats.hotspotMutationCount());
    }
}