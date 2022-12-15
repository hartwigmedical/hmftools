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
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.orange.TestOrangeReportFactory;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.linx.TestLinxInterpretationFactory;
import com.hartwig.hmftools.orange.algo.purple.ImmutablePurityPloidyFit;
import com.hartwig.hmftools.orange.algo.purple.PurityPloidyFit;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurpleVariant;
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
                        .allGermlineDeletions(Lists.newArrayList())
                        .reportableGermlineDeletions(Lists.newArrayList())
                        .build())
                .linx(TestLinxInterpretationFactory.builder()
                        .allGermlineDisruptions(Lists.newArrayList())
                        .reportableGermlineDisruptions(Lists.newArrayList())
                        .build())
                .build();

        OrangeReport converted = GermlineConversion.convertGermlineToSomatic(report);

        assertNull(converted.germlineMVLHPerGene());

        assertNull(converted.purple().germlineDrivers());
        assertNull(converted.purple().allGermlineVariants());
        assertNull(converted.purple().reportableGermlineVariants());
        assertNull(converted.purple().additionalSuspectGermlineVariants());
        assertNull(converted.purple().allGermlineDeletions());
        assertNull(converted.purple().reportableGermlineDeletions());

        assertNull(converted.linx().allGermlineDisruptions());
        assertNull(converted.linx().reportableGermlineDisruptions());
    }

    @Test
    public void canConvertPurple() {
        PurpleVariant somaticVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant reportableSomaticVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant suspectSomaticVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant germlineVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant reportableGermlineVariant = TestPurpleVariantFactory.builder().build();
        PurpleVariant suspectGermlineVariant = TestPurpleVariantFactory.builder().build();

        DriverCatalog somaticDriver = DriverCatalogTestFactory.builder().driver(DriverType.AMP).build();
        DriverCatalog germlineMutationDriver = DriverCatalogTestFactory.builder().driver(DriverType.GERMLINE_MUTATION).build();
        DriverCatalog germlineDisruptionDriver = DriverCatalogTestFactory.builder().driver(DriverType.GERMLINE_DISRUPTION).build();

        PurpleInterpretedData purple = TestPurpleInterpretationFactory.builder()
                .fit(createWithGermlineAberration())
                .addSomaticDrivers(somaticDriver)
                .addGermlineDrivers(germlineMutationDriver, germlineDisruptionDriver)
                .addAllSomaticVariants(somaticVariant)
                .addReportableSomaticVariants(reportableSomaticVariant)
                .addAdditionalSuspectSomaticVariants(suspectSomaticVariant)
                .addAllGermlineVariants(germlineVariant)
                .addReportableGermlineVariants(reportableGermlineVariant)
                .addAdditionalSuspectGermlineVariants(suspectGermlineVariant)
                .build();

        PurpleInterpretedData converted = GermlineConversion.convertPurpleGermline(purple);

        assertTrue(converted.fit().qc().germlineAberrations().isEmpty());

        assertNull(converted.germlineDrivers());
        assertNull(converted.allGermlineVariants());
        assertNull(converted.reportableGermlineVariants());
        assertNull(converted.additionalSuspectGermlineVariants());

        assertEquals(2, converted.somaticDrivers().size());
        assertNotNull(findByDriverType(converted.somaticDrivers(), DriverType.AMP));
        assertNotNull(findByDriverType(converted.somaticDrivers(), DriverType.MUTATION));

        assertEquals(2, converted.allSomaticVariants().size());
        assertTrue(converted.allSomaticVariants().contains(somaticVariant));
        assertTrue(converted.allSomaticVariants().contains(germlineVariant));

        assertEquals(2, converted.reportableSomaticVariants().size());
        assertTrue(converted.reportableSomaticVariants().contains(reportableSomaticVariant));
        assertTrue(converted.reportableSomaticVariants().contains(reportableGermlineVariant));

        assertEquals(1, converted.additionalSuspectSomaticVariants().size());
        assertTrue(converted.additionalSuspectSomaticVariants().contains(suspectSomaticVariant));
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
    public void canMergeDrivers() {
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

}