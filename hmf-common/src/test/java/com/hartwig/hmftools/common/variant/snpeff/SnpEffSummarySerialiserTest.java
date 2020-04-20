package com.hartwig.hmftools.common.variant.snpeff;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Random;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.junit.Test;

public class SnpEffSummarySerialiserTest {

    @Test
    public void testSerialise() {
        SnpEffSummary random = createRandom(new Random()).build();
        List<String> worst = SnpEffSummarySerialiser.worstDetails(random);
        List<String> canonical = SnpEffSummarySerialiser.canonicalDetails(random);
        SnpEffSummary recreated = SnpEffSummarySerialiser.fromDetails(worst, canonical);

        assertEquals(random, recreated);
    }

    @Test
    public void testEffect() {
        SnpEffSummary random = createRandom(new Random()).worstEffect("multiple effects; with spaces").build();
        List<String> worst = SnpEffSummarySerialiser.worstDetails(random);
        assertEquals("multiple_effects&with_spaces", worst.get(2));
    }

    static ImmutableSnpEffSummary.Builder createRandom(Random random) {
        return ImmutableSnpEffSummary.builder()
                .genesAffected(random.nextInt())
                .worstGene(Integer.toString(random.nextInt()))
                .worstEffect(Integer.toString(random.nextInt()))
                .worstCodingEffect(CodingEffect.values()[random.nextInt(CodingEffect.values().length)])
                .worstTranscript(Integer.toString(random.nextInt()))
                .canonicalGene(Integer.toString(random.nextInt()))
                .canonicalEffect(Integer.toString(random.nextInt()))
                .canonicalTranscript(Integer.toString(random.nextInt()))
                .canonicalCodingEffect(CodingEffect.values()[random.nextInt(CodingEffect.values().length)])
                .canonicalHgvsCodingImpact(Integer.toString(random.nextInt()))
                .canonicalHgvsProteinImpact(Integer.toString(random.nextInt()));
    }
}
