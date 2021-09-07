package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DriverCatalogFileTest {

    private static final String DRIVER_CATALOG_TSV = Resources.getResource("drivercatalog/sample.driver.catalog.tsv").getPath();
    private static final double EPSILON = 1.0e-10;

    @Test
    public void canReadDriverCatalogFileWithoutDnds() throws IOException {
        List<DriverCatalog> driverCatalogList = DriverCatalogFile.read(DRIVER_CATALOG_TSV);
        assertEquals(3, driverCatalogList.size());

        DriverCatalog catalog1 = driverCatalogList.get(0);
        assertEquals("8", catalog1.chromosome());
        assertEquals("p22", catalog1.chromosomeBand());
        assertEquals("SGCZ", catalog1.gene());
        assertEquals(DriverCategory.TSG, catalog1.category());
        assertEquals(DriverType.HOM_DISRUPTION, catalog1.driver());
        assertEquals(LikelihoodMethod.DNDS, catalog1.likelihoodMethod());
        assertEquals(1.0000, catalog1.driverLikelihood(), EPSILON);
        assertEquals(0, catalog1.missense(), EPSILON);
        assertEquals(0, catalog1.nonsense(), EPSILON);
        assertEquals(0, catalog1.splice(), EPSILON);
        assertEquals(0, catalog1.inframe(), EPSILON);
        assertEquals(0, catalog1.frameshift(), EPSILON);
        assertFalse(catalog1.biallelic());
        assertEquals(3.7, catalog1.minCopyNumber(), EPSILON);
        assertEquals(1.3, catalog1.maxCopyNumber(), EPSILON);

        DriverCatalog catalog2 = driverCatalogList.get(1);
        assertEquals("X", catalog2.chromosome());
        assertEquals("q21.33", catalog2.chromosomeBand());
        assertEquals("DIAPH2", catalog2.gene());
        assertEquals(DriverCategory.ONCO, catalog2.category());
        assertEquals(DriverType.MUTATION, catalog2.driver());
        assertEquals(LikelihoodMethod.BIALLELIC, catalog2.likelihoodMethod());
        assertEquals(1.0000, catalog2.driverLikelihood(), EPSILON);
        assertEquals(0, catalog2.missense(), EPSILON);
        assertEquals(0, catalog2.nonsense(), EPSILON);
        assertEquals(0, catalog2.splice(), EPSILON);
        assertEquals(0, catalog2.inframe(), EPSILON);
        assertEquals(0, catalog2.frameshift(), EPSILON);
        assertTrue(catalog2.biallelic());
        assertEquals(2.4, catalog2.minCopyNumber(), EPSILON);
        assertEquals(5.3, catalog2.maxCopyNumber(), EPSILON);

        DriverCatalog catalog3 = driverCatalogList.get(2);
        assertEquals("9", catalog3.chromosome());
        assertEquals("p23-p24.1", catalog3.chromosomeBand());
        assertEquals("PTPRD", catalog3.gene());
        assertEquals(DriverCategory.TSG, catalog3.category());
        assertEquals(DriverType.HOM_DISRUPTION, catalog3.driver());
        assertEquals(LikelihoodMethod.DEL, catalog3.likelihoodMethod());
        assertEquals(1.0000, catalog3.driverLikelihood(), EPSILON);
        assertEquals(0, catalog3.missense(), EPSILON);
        assertEquals(0, catalog3.nonsense(), EPSILON);
        assertEquals(0, catalog3.splice(), EPSILON);
        assertEquals(0, catalog3.inframe(), EPSILON);
        assertEquals(0, catalog3.frameshift(), EPSILON);
        assertFalse(catalog3.biallelic());
        assertEquals(0.8, catalog3.minCopyNumber(), EPSILON);
        assertEquals(2.4, catalog3.maxCopyNumber(), EPSILON);
    }

    @Test
    public void canSerializeDriverCatalogToStringAndBack() {
        final Random random = new Random();
        final DriverCatalog input = createRandomScore(random);

        final List<String> lines = DriverCatalogFile.toLines(Lists.newArrayList(input));
        assertEquals(2, lines.size());
        assertTrue(lines.get(0).startsWith("chromosome"));

        assertEquals(input, DriverCatalogFile.fromLines(lines).get(0));
    }

    @NotNull
    private static DriverCatalog createRandomScore(@NotNull Random random) {
        return ImmutableDriverCatalog.builder()
                .chromosome("" + random.nextLong())
                .chromosomeBand("" + random.nextLong())
                .gene("" + random.nextLong())
                .driver(DriverType.values()[random.nextInt(DriverType.values().length)])
                .category(DriverCategory.values()[random.nextInt(DriverCategory.values().length)])
                .likelihoodMethod(LikelihoodMethod.values()[random.nextInt(LikelihoodMethod.values().length)])
                .driverLikelihood(nextDouble(random))
                .missense(random.nextLong())
                .nonsense(random.nextLong())
                .splice(random.nextLong())
                .inframe(random.nextLong())
                .frameshift(random.nextLong())
                .biallelic(random.nextBoolean())
                .minCopyNumber(nextDouble(random))
                .maxCopyNumber(nextDouble(random))
                .build();
    }

    private static double nextDouble(@NotNull final Random random) {
        return Double.parseDouble(DriverCatalogFile.FORMAT.format(random.nextDouble()));
    }
}
