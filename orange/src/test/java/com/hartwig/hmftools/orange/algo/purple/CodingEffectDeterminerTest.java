package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantEffect;

import org.junit.Test;

public class CodingEffectDeterminerTest
{
    @Test
    public void testEmptyEffectListReturnsNoneCodingEffect()
    {
        List<VariantEffect> emptyVariantEffectList = List.of();

        CodingEffect actualCodingEffect = CodingEffectDeterminer.determineCodingEffect(emptyVariantEffectList);

        assertEquals(CodingEffect.NONE, actualCodingEffect);
    }

    @Test
    public void testNonsenseOrFrameshiftDetermination()
    {
        List<VariantEffect> variantEffects = List.of(VariantEffect.START_LOST);

        CodingEffect actualCodingEffect = CodingEffectDeterminer.determineCodingEffect(variantEffects);

        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, actualCodingEffect);
    }

    @Test
    public void testSpliceCodingEffectDetermination()
    {
        List<VariantEffect> variantEffects = List.of(VariantEffect.SPLICE_DONOR);

        CodingEffect actualCodingEffect = CodingEffectDeterminer.determineCodingEffect(variantEffects);

        assertEquals(CodingEffect.SPLICE, actualCodingEffect);
    }

    @Test
    public void testMissenseEffectDetermination()
    {
        List<VariantEffect> variantEffects = List.of(VariantEffect.MISSENSE);

        CodingEffect actualCodingEffect = CodingEffectDeterminer.determineCodingEffect(variantEffects);

        assertEquals(CodingEffect.MISSENSE, actualCodingEffect);
    }

    @Test
    public void testSynonymousEffectDetermination()
    {
        List<VariantEffect> variantEffects = List.of(VariantEffect.SYNONYMOUS);

        CodingEffect actualCodingEffect = CodingEffectDeterminer.determineCodingEffect(variantEffects);

        assertEquals(CodingEffect.SYNONYMOUS, actualCodingEffect);
    }

    @Test
    public void testOtherEffectDetermination()
    {
        List<VariantEffect> variantEffects = List.of(VariantEffect.INTRONIC);

        CodingEffect actualCodingEffect = CodingEffectDeterminer.determineCodingEffect(variantEffects);

        assertEquals(CodingEffect.NONE, actualCodingEffect);
    }

    @Test
    public void nonsenseOrFrameshiftTakesPrecedenceOverSplice()
    {
        List<VariantEffect> variantEffects = List.of(VariantEffect.START_LOST, VariantEffect.SPLICE_DONOR);

        CodingEffect actualCodingEffect = CodingEffectDeterminer.determineCodingEffect(variantEffects);

        assertEquals(CodingEffect.NONSENSE_OR_FRAMESHIFT, actualCodingEffect);
    }

    @Test
    public void spliceTakesPrecedenceOverMissense()
    {
        List<VariantEffect> variantEffects = List.of(VariantEffect.MISSENSE, VariantEffect.SPLICE_DONOR);

        CodingEffect actualCodingEffect = CodingEffectDeterminer.determineCodingEffect(variantEffects);

        assertEquals(CodingEffect.SPLICE, actualCodingEffect);
    }

    @Test
    public void missenseTakesPrecedenceOverSynonymous()
    {
        List<VariantEffect> variantEffects = List.of(VariantEffect.MISSENSE, VariantEffect.SYNONYMOUS);

        CodingEffect actualCodingEffect = CodingEffectDeterminer.determineCodingEffect(variantEffects);

        assertEquals(CodingEffect.MISSENSE, actualCodingEffect);
    }

    @Test
    public void synonymousTakesPrecedenceOverNone()
    {
        List<VariantEffect> variantEffects = List.of(VariantEffect.INTRONIC, VariantEffect.SYNONYMOUS);

        CodingEffect actualCodingEffect = CodingEffectDeterminer.determineCodingEffect(variantEffects);

        assertEquals(CodingEffect.SYNONYMOUS, actualCodingEffect);
    }
}
