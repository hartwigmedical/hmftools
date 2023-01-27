package com.hartwig.hmftools.orange.algo.util;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogTestFactory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.orange.TestOrangeReportFactory;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpretedData;
import com.hartwig.hmftools.orange.algo.linx.TestLinxInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.ImmutablePurityPloidyFit;
import com.hartwig.hmftools.orange.algo.purple.PurityPloidyFit;
import com.hartwig.hmftools.orange.algo.purple.PurpleGainLoss;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleGainLossFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.TestPurpleVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class GermlineConversionTest {

    private static final double EPSILON = 1.0E-10;

    @Test
    public void wipesAllGermlineUponConversion() {
        OrangeReport report = TestOrangeReportFactory.builder()
                .germlineMVLHPerGene(Maps.newHashMap())
                .purple(TestPurpleInterpretationFactory.builder()
                        .germlineDrivers(Lists.newArrayList())
                        .allGermlineVariants(Lists.newArrayList())
                        .reportableGermlineVariants(Lists.newArrayList())
                        .additionalSuspectGermlineVariants(Lists.newArrayList())
                        .allGermlineGainsLosses(Lists.newArrayList())
                        .reportableGermlineGainsLosses(Lists.newArrayList())
                        .build())
                .linx(TestLinxInterpretationFactory.builder()
                        .allGermlineStructuralVariants(Lists.newArrayList())
                        .allGermlineBreakends(Lists.newArrayList())
                        .reportableGermlineBreakends(Lists.newArrayList())
                        .build())
                .build();

        OrangeReport converted = GermlineConversion.convertGermlineToSomatic(report);

        assertNull(converted.germlineMVLHPerGene());

        assertNull(converted.purple().germlineDrivers());
        assertNull(converted.purple().allGermlineVariants());
        assertNull(converted.purple().reportableGermlineVariants());
        assertNull(converted.purple().additionalSuspectGermlineVariants());
        assertNull(converted.purple().allGermlineGainsLosses());
        assertNull(converted.purple().reportableGermlineGainsLosses());

        assertNull(converted.linx().allGermlineStructuralVariants());
        assertNull(converted.linx().allGermlineBreakends());
        assertNull(converted.linx().reportableGermlineBreakends());
    }

    @Test
    public void canConvertMinimalPurpleData() {
        OrangeReport minimal = TestOrangeReportFactory.createMinimalTestReport();

        assertNotNull(GermlineConversion.convertPurpleGermline(true, minimal.purple()));
    }

    @Test
    public void canConvertPurple() {
        PurpleVariant somaticVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant reportableSomaticVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant suspectSomaticVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant germlineVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant reportableGermlineVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant suspectGermlineVariant = TestPurpleVariantFactory.builder().build();

        PurpleGainLoss somaticGainLoss = TestPurpleGainLossFactory.builder().build();
        PurpleGainLoss reportableSomaticGainLoss = TestPurpleGainLossFactory.builder().build();
        PurpleGainLoss germlineGainLoss = TestPurpleGainLossFactory.builder().build();
        PurpleGainLoss reportableGermlineGainLoss = TestPurpleGainLossFactory.builder().build();

        DriverCatalog somaticDriver = DriverCatalogTestFactory.builder().driver(DriverType.AMP).build();
        DriverCatalog germlineMutationDriver = DriverCatalogTestFactory.builder().driver(DriverType.GERMLINE_MUTATION).build();
        DriverCatalog germlineDisruptionDriver = DriverCatalogTestFactory.builder().driver(DriverType.GERMLINE_DISRUPTION).build();

        PurpleInterpretedData purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineMutationDriver, germlineDisruptionDriver)
                .addAllSomaticVariants(somaticVariant, suspectSomaticVariant, reportableSomaticVariant)
                .addReportableSomaticVariants(reportableSomaticVariant)
                .addAdditionalSuspectSomaticVariants(suspectSomaticVariant)
                .addAllGermlineVariants(germlineVariant, suspectGermlineVariant, reportableGermlineVariant)
                .addReportableGermlineVariants(reportableGermlineVariant)
                .addAdditionalSuspectGermlineVariants(suspectGermlineVariant)
                .addAllSomaticGainsLosses(somaticGainLoss, reportableSomaticGainLoss)
                .addReportableSomaticGainsLosses(reportableSomaticGainLoss)
                .addAllGermlineGainsLosses(germlineGainLoss, reportableGermlineGainLoss)
                .addReportableGermlineGainsLosses(reportableGermlineGainLoss)
                .build();

        PurpleInterpretedData converted = GermlineConversion.convertPurpleGermline(true, purple);

        assertTrue(converted.fit().qc().germlineAberrations().isEmpty());

        assertEquals(3, converted.somaticDrivers().size());
        assertNotNull(findByDriverType(converted.somaticDrivers(), DriverType.AMP));
        assertNotNull(findByDriverType(converted.somaticDrivers(), DriverType.MUTATION));
        assertNotNull(findByDriverType(converted.somaticDrivers(), DriverType.DISRUPTION));

        assertEquals(4, converted.allSomaticVariants().size());
        assertTrue(converted.allSomaticVariants().contains(reportableGermlineVariant));

        assertEquals(2, converted.reportableSomaticVariants().size());
        assertTrue(converted.reportableSomaticVariants().contains(reportableSomaticVariant));
        assertTrue(converted.reportableSomaticVariants().contains(reportableGermlineVariant));

        assertEquals(1, converted.additionalSuspectSomaticVariants().size());
        assertTrue(converted.additionalSuspectSomaticVariants().contains(suspectSomaticVariant));

        assertEquals(3, converted.allSomaticGainsLosses().size());
        assertTrue(converted.allSomaticGainsLosses().contains(reportableGermlineGainLoss));

        assertEquals(2, converted.reportableSomaticGainsLosses().size());
        assertTrue(converted.reportableSomaticGainsLosses().contains(reportableSomaticGainLoss));
        assertTrue(converted.reportableSomaticGainsLosses().contains(reportableGermlineGainLoss));

        PurpleInterpretedData unreliableConverted = GermlineConversion.convertPurpleGermline(false, purple);
        assertEquals(1, unreliableConverted.somaticDrivers().size());
        assertNotNull(findByDriverType(unreliableConverted.somaticDrivers(), DriverType.AMP));

        assertEquals(1, unreliableConverted.reportableSomaticVariants().size());
        assertTrue(unreliableConverted.reportableSomaticVariants().contains(reportableSomaticVariant));

        assertEquals(1, unreliableConverted.reportableSomaticGainsLosses().size());
        assertTrue(unreliableConverted.reportableSomaticGainsLosses().contains(reportableSomaticGainLoss));
    }

    @NotNull
    private static PurityPloidyFit createWithGermlineAberration() {
        PurityPloidyFit base = TestPurpleInterpretationFactory.createMinimalTestFitData();
        return ImmutablePurityPloidyFit.builder()
                .from(base)
                .qc(ImmutablePurpleQC.builder().from(base.qc()).addGermlineAberrations(GermlineAberration.MOSAIC_X).build())
                .build();
    }

    @Nullable
    private static DriverCatalog findByDriverType(@NotNull List<DriverCatalog> drivers, @NotNull DriverType driverTypeToFind) {
        for (DriverCatalog driver : drivers) {
            if (driver.driver() == driverTypeToFind) {
                return driver;
            }
        }

        return null;
    }

    @Test
    public void canMergeMutationDrivers() {
        DriverCatalog somaticDriver1 = DriverCatalogTestFactory.builder()
                .driver(DriverType.MUTATION)
                .gene("gene 1")
                .transcript("transcript 1")
                .driverLikelihood(0.5)
                .likelihoodMethod(LikelihoodMethod.DNDS)
                .missense(1)
                .nonsense(2)
                .splice(3)
                .inframe(4)
                .frameshift(5)
                .biallelic(false)
                .build();

        DriverCatalog somaticDriver2 =
                DriverCatalogTestFactory.builder().driver(DriverType.MUTATION).gene("gene 1").transcript("transcript 2").build();

        DriverCatalog germlineDriver1 = DriverCatalogTestFactory.builder()
                .driver(DriverType.GERMLINE_MUTATION)
                .gene("gene 1")
                .transcript("transcript 1")
                .driverLikelihood(0.8)
                .likelihoodMethod(LikelihoodMethod.GERMLINE)
                .missense(1)
                .nonsense(2)
                .splice(3)
                .inframe(4)
                .frameshift(5)
                .biallelic(true)
                .build();

        DriverCatalog germlineDriver2 =
                DriverCatalogTestFactory.builder().driver(DriverType.GERMLINE_MUTATION).gene("gene 2").transcript("transcript 1").build();

        List<DriverCatalog> merged = GermlineConversion.mergeGermlineDriversIntoSomatic(Lists.newArrayList(somaticDriver1, somaticDriver2),
                Lists.newArrayList(germlineDriver1, germlineDriver2));

        assertEquals(3, merged.size());
        DriverCatalog driver1 = findByGeneTranscript(merged, "gene 1", "transcript 1");
        assertNotNull(driver1);
        assertEquals(DriverType.MUTATION, driver1.driver());
        assertEquals(LikelihoodMethod.HOTSPOT, driver1.likelihoodMethod());
        assertEquals(0.8, driver1.driverLikelihood(), EPSILON);
        assertEquals(2, driver1.missense());
        assertEquals(4, driver1.nonsense());
        assertEquals(6, driver1.splice());
        assertEquals(8, driver1.inframe());
        assertEquals(10, driver1.frameshift());
        assertTrue(driver1.biallelic());

        DriverCatalog driver2 = findByGeneTranscript(merged, "gene 1", "transcript 2");
        assertEquals(somaticDriver2, driver2);

        DriverCatalog driver3 = findByGeneTranscript(merged, "gene 2", "transcript 1");
        assertEquals(DriverType.MUTATION, driver3.driver());
        assertEquals(LikelihoodMethod.HOTSPOT, driver3.likelihoodMethod());
    }

    @Test
    public void canConvertAllTypesOfGermlineDrivers() {
        DriverCatalog somaticDriver =
                DriverCatalogTestFactory.builder().gene("gene 1").transcript("transcript 1").driver(DriverType.DEL).build();

        List<DriverCatalog> mergedNoGermline = GermlineConversion.mergeGermlineDriversIntoSomatic(Lists.newArrayList(somaticDriver), null);
        assertEquals(1, mergedNoGermline.size());
        assertTrue(mergedNoGermline.contains(somaticDriver));

        DriverCatalog germlineDriver1 =
                DriverCatalogTestFactory.builder().gene("gene 1").transcript("transcript 1").driver(DriverType.GERMLINE_DELETION).build();

        DriverCatalog germlineDriver2 =
                DriverCatalogTestFactory.builder().gene("gene 2").transcript("transcript 2").driver(DriverType.GERMLINE_DELETION).build();

        DriverCatalog germlineDriver3 =
                DriverCatalogTestFactory.builder().gene("gene 3").transcript("transcript 3").driver(DriverType.GERMLINE_MUTATION).build();

        DriverCatalog germlineDriver4 = DriverCatalogTestFactory.builder()
                .gene("gene 4")
                .transcript("transcript 4")
                .driver(DriverType.GERMLINE_DISRUPTION)
                .driverLikelihood(1D)
                .build();

        List<DriverCatalog> merged = GermlineConversion.mergeGermlineDriversIntoSomatic(Lists.newArrayList(somaticDriver),
                Lists.newArrayList(germlineDriver1, germlineDriver2, germlineDriver3, germlineDriver4));

        assertEquals(4, merged.size());
        assertTrue(mergedNoGermline.contains(somaticDriver));

        DriverCatalog germlineDeletionDriver = findByGeneTranscript(merged, germlineDriver2.gene(), germlineDriver2.transcript());
        assertEquals(DriverType.DEL, germlineDeletionDriver.driver());

        DriverCatalog germlineMutationDriver = findByGeneTranscript(merged, germlineDriver3.gene(), germlineDriver3.transcript());
        assertEquals(DriverType.MUTATION, germlineMutationDriver.driver());

        DriverCatalog germlineDisruptionDriver = findByGeneTranscript(merged, germlineDriver4.gene(), germlineDriver4.transcript());
        assertEquals(DriverType.DISRUPTION, germlineDisruptionDriver.driver());
        assertEquals(0D, germlineDisruptionDriver.driverLikelihood(), EPSILON);
    }

    @Nullable
    private static DriverCatalog findByGeneTranscript(@NotNull List<DriverCatalog> drivers, @NotNull String geneToFind,
            @NotNull String transcriptToFind) {
        for (DriverCatalog driver : drivers) {
            if (driver.gene().equals(geneToFind) && driver.transcript().equals(transcriptToFind)) {
                return driver;
            }
        }
        return null;
    }

    @Test
    public void canConvertLinx() {
        LinxSvAnnotation somaticStructuralVariant1 = LinxTestFactory.svAnnotationBuilder().svId(1).clusterId(5).build();
        LinxSvAnnotation somaticStructuralVariant2 = LinxTestFactory.svAnnotationBuilder().svId(2).clusterId(6).build();
        LinxBreakend somaticBreakend = LinxTestFactory.breakendBuilder().id(8).svId(1).build();
        LinxBreakend reportableSomaticBreakend = LinxTestFactory.breakendBuilder().id(9).svId(2).build();

        LinxSvAnnotation germlineStructuralVariant1 = LinxTestFactory.svAnnotationBuilder().svId(1).clusterId(5).build();
        LinxSvAnnotation germlineStructuralVariant2 = LinxTestFactory.svAnnotationBuilder().svId(2).clusterId(6).build();
        LinxBreakend germlineBreakend = LinxTestFactory.breakendBuilder().id(8).svId(1).build();
        LinxBreakend reportableGermlineBreakend = LinxTestFactory.breakendBuilder().id(9).svId(2).build();

        LinxInterpretedData linx = TestLinxInterpretationFactory.builder()
                .addAllSomaticStructuralVariants(somaticStructuralVariant1, somaticStructuralVariant2)
                .addAllSomaticBreakends(somaticBreakend, reportableSomaticBreakend)
                .addReportableSomaticBreakends(reportableSomaticBreakend)
                .addAllGermlineStructuralVariants(germlineStructuralVariant1, germlineStructuralVariant2)
                .addAllGermlineBreakends(germlineBreakend, reportableGermlineBreakend)
                .addReportableGermlineBreakends(reportableGermlineBreakend)
                .build();

        LinxInterpretedData converted = GermlineConversion.convertLinxGermline(true, linx);
        assertEquals(4, converted.allSomaticStructuralVariants().size());
        assertEquals(4, GermlineConversion.findMaxSvId(converted.allSomaticStructuralVariants()));
        assertEquals(8, GermlineConversion.findMaxClusterId(converted.allSomaticStructuralVariants()));

        assertEquals(3, converted.allSomaticBreakends().size());
        assertEquals(11, GermlineConversion.findMaxBreakendId(converted.allSomaticBreakends()));

        assertEquals(2, converted.reportableSomaticBreakends().size());
        assertEquals(11, GermlineConversion.findMaxBreakendId(converted.reportableSomaticBreakends()));
    }
}