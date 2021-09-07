package com.hartwig.hmftools.common.variant.snpeff;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Random;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class VariantImpactSerialiserTest
{
    @Test
    public void testSerialise()
    {
        Random random = new Random();

        VariantImpact variantImpact = new VariantImpact(
                random.nextInt(), "", Integer.toString(random.nextInt()), Integer.toString(random.nextInt()), Integer.toString(random.nextInt()),
                CodingEffect.values()[random.nextInt(CodingEffect.values().length)], Integer.toString(random.nextInt()),
                Integer.toString(random.nextInt()), Integer.toString(random.nextInt()), Integer.toString(random.nextInt()),
                Integer.toString(random.nextInt()), CodingEffect.values()[random.nextInt(CodingEffect.values().length)]);

        List<String> worst = VariantImpactSerialiser.worstDetails(variantImpact);
        List<String> canonical = VariantImpactSerialiser.canonicalDetails(variantImpact);
        VariantImpact recreated = VariantImpactSerialiser.fromVcfAnnotation(worst, canonical);

        // TODO - add equality method
        // assertEquals(variantImpact, recreated);
    }

    @Test
    public void testEffect()
    {
        Random random = new Random();

        VariantImpact variantImpact = new VariantImpact(
                random.nextInt(), "", Integer.toString(random.nextInt()), Integer.toString(random.nextInt()), Integer.toString(random.nextInt()),
                CodingEffect.values()[random.nextInt(CodingEffect.values().length)], Integer.toString(random.nextInt()),
                Integer.toString(random.nextInt()), Integer.toString(random.nextInt()), "multiple effects; with spaces",
                Integer.toString(random.nextInt()), CodingEffect.values()[random.nextInt(CodingEffect.values().length)]);

        List<String> worst = VariantImpactSerialiser.worstDetails(variantImpact);
        assertEquals("multiple_effects&with_spaces", worst.get(2));
    }
}
