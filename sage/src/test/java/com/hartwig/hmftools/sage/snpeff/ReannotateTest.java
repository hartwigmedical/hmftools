package com.hartwig.hmftools.sage.snpeff;

import static java.util.Collections.singletonList;

import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_DELETION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INFRAME_INSERTION;
import static com.hartwig.hmftools.common.variant.VariantConsequence.sufferConsequences;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.variant.CodingEffect;

import org.junit.Test;

public class ReannotateTest {

    @Test
    public void testInframeIndelStringsProduceExpectedCodingEffect() {
        assertEquals(INFRAME_INSERTION, sufferConsequences(singletonList(Reannotate.INFRAME_INSERTION)).get(0));
        assertEquals(INFRAME_DELETION, sufferConsequences(singletonList(Reannotate.INFRAME_DELETION)).get(0));

        assertEquals(CodingEffect.MISSENSE, CodingEffect.effect("Gene", singletonList(INFRAME_INSERTION)));
        assertEquals(CodingEffect.MISSENSE, CodingEffect.effect("Gene", singletonList(INFRAME_DELETION)));
    }

}
