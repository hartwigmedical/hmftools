package com.hartwig.hmftools.common.drivercatalog.dnds;

import static org.junit.Assert.assertEquals;

import java.util.Random;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class DndsVariantFileTest {

    @Test
    public void testSerialization() {
        DndsVariant victim = random(new Random());
        assertEquals(victim, DndsVariantFile.fromString(DndsVariantFile.toString(victim)));
    }

    @NotNull
    private static DndsVariant random(@NotNull final Random random) {
        return ImmutableDndsVariant.builder()
                .sampleId("" + random.nextInt())
                .chromosome("" + random.nextInt())
                .position(random.nextInt())
                .ref("" + random.nextInt())
                .alt("" + random.nextInt())
                .gene("" + random.nextInt())
                .worstCodingEffect(CodingEffect.values()[random.nextInt(CodingEffect.values().length - 1)])
                .canonicalCodingEffect(CodingEffect.values()[random.nextInt(CodingEffect.values().length - 1)])
                .biallelic(random.nextBoolean())
                .hotspot(random.nextBoolean())
                .repeatCount(random.nextInt())
                .build();

    }
}
