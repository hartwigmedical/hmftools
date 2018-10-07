package com.hartwig.hmftools.patientreporter.variants;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.ImmutableEnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantTestBuilderFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class SomaticVariantAnalyzerTest {

    private static final String PASS_FILTER = "PASS";

    private static final CodingEffect RIGHT_EFFECT = CodingEffect.MISSENSE;
    private static final CodingEffect WRONG_EFFECT = CodingEffect.SYNONYMOUS;
    private static final String RIGHT_GENE = "RIGHT";
    private static final String WRONG_GENE = "WRONG";

    @Test
    public void realCaseWorks() {
        final List<EnrichedSomaticVariant> variants =
                Lists.newArrayList(builder().gene(RIGHT_GENE).canonicalCodingEffect(RIGHT_EFFECT).worstCodingEffect(RIGHT_EFFECT).build(),
                        builder().gene(RIGHT_GENE).canonicalCodingEffect(WRONG_EFFECT).worstCodingEffect(WRONG_EFFECT).build(),
                        builder().gene(WRONG_GENE).canonicalCodingEffect(RIGHT_EFFECT).worstCodingEffect(RIGHT_EFFECT).build(),
                        builder().gene(WRONG_GENE).canonicalCodingEffect(WRONG_EFFECT).worstCodingEffect(WRONG_EFFECT).build());

        final SomaticVariantAnalysis analysis = SomaticVariantAnalyzer.run(variants, Sets.newHashSet(RIGHT_GENE));

        assertEquals(2, analysis.tumorMutationalLoad());
        assertEquals(1, analysis.variantsToReport().size());
    }

    @NotNull
    private static ImmutableEnrichedSomaticVariant.Builder builder() {
        return SomaticVariantTestBuilderFactory.createEnriched().filter(PASS_FILTER);
    }
}
