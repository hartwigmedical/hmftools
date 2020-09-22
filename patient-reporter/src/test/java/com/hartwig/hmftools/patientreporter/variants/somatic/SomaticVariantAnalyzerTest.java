package com.hartwig.hmftools.patientreporter.variants.somatic;

import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestDriverCatalogBuilder;
import static com.hartwig.hmftools.patientreporter.PatientReporterTestFactory.createTestSomaticVariantBuilder;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.junit.Test;

public class SomaticVariantAnalyzerTest {

    @Test
    public void canAnalyzeSomaticVariants() {
        double driverLikelihood = 0.8;

        SomaticVariant variant = createTestSomaticVariantBuilder().gene("GENE").reported(true).build();
        DriverCatalog driver = createTestDriverCatalogBuilder().gene("GENE").driverLikelihood(driverLikelihood).build();

        List<DriverSomaticVariant> drivers = SomaticVariantAnalyzer.run(Lists.newArrayList(variant), Lists.newArrayList(driver));
        assertEquals(1, drivers.size());
        assertEquals(variant, drivers.get(0).variant());
        assertEquals(driverLikelihood, drivers.get(0).driverLikelihood(), 1.0E-10);
    }

    @Test (expected = IllegalStateException.class)
    public void crashOnMissingDriverCatalogEntry() {
        SomaticVariant variant = createTestSomaticVariantBuilder().gene("GENE").reported(true).build();

        SomaticVariantAnalyzer.run(Lists.newArrayList(variant), Lists.newArrayList());
    }

}