package com.hartwig.hmftools.common.drivercatalog;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DriverCatalogFileTest {

    @Test
    public void testInputAndOutput() {
        final Random random = new Random();
        final DriverCatalog input = createRandomScore(random);

        final List<String> lines = DriverCatalogFile.toLines(Lists.newArrayList(input));
        assertEquals(2, lines.size());
        assertTrue(lines.get(0).startsWith(DriverCatalogFile.HEADER_PREFIX));

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
                .dndsLikelihood(nextDouble(random))
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
        return Double.valueOf(DriverCatalogFile.FORMAT.format(random.nextDouble()));
    }
}
