package com.hartwig.hmftools.vian;

import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

public class RefVariantData
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public final VariantType Type;
    public final String Gene;

    public final String CanonicalEffect;
    public final CodingEffect CanonicalCodingEffect;
    public final String CanonicalHgvsCodingImpact;
    public final String CanonicalHgvsProteinImpact;

    public final String WorstEffect;
    public final CodingEffect WorstCodingEffect;
    public final String WorstEffectTranscript;
    public final int GenesAffected;

    public final String Microhomology;
    public final String RepeatSequence;
    public final boolean PhasedInframeIndel;

    // canonicalHgvsProteinImpact, , repeatCount, localPhaseSet, localRealignmentSet, phasedInframeIndel, reported

    public RefVariantData(
            final String chromosome, final int position, final String ref, final String alt, final VariantType type,
            final String gene, final String canonicalEffect, final CodingEffect canonicalCodingEffect,
            final String worstEffect, final CodingEffect worstCodingEffect, final String worstEffectTranscript, final int genesAffected,
            final String canonicalHgvsCodingImpact, final String canonicalHgvsProteinImpact,
            final String microhomology, final String repeatSequence, boolean phasedInframeIndel)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
        Gene = gene;

        CanonicalEffect = canonicalEffect;
        CanonicalCodingEffect = canonicalCodingEffect;
        WorstEffect = worstEffect;
        WorstCodingEffect = worstCodingEffect;
        WorstEffectTranscript = worstEffectTranscript;
        GenesAffected = genesAffected;

        CanonicalHgvsCodingImpact = canonicalHgvsCodingImpact;
        CanonicalHgvsProteinImpact = canonicalHgvsProteinImpact;
        Microhomology = microhomology;
        RepeatSequence = repeatSequence;
        PhasedInframeIndel = phasedInframeIndel;
    }

    public static RefVariantData fromSomatic(final SomaticVariant variant)
    {
        return new RefVariantData(
                variant.chromosome(), (int)variant.position(), variant.ref(), variant.alt(), variant.type(), variant.gene(),
                variant.canonicalEffect(), variant.canonicalCodingEffect(),
                variant.worstEffect(), variant.worstCodingEffect(), variant.worstEffectTranscript(),
                variant.genesAffected(),  variant.canonicalHgvsCodingImpact(), variant.canonicalHgvsProteinImpact(),
                variant.microhomology(), variant.repeatSequence(), variant.phasedInframeIndelIdentifier() > 0);
    }

    public String toString()
    {
        return String.format("pos(%s:%d) variant(%s: %s>%s) canon(%s: %s) worst(%s: %s)",
                Chromosome, Position, Type, Ref, Alt, CanonicalCodingEffect, CanonicalEffect,
                WorstCodingEffect, WorstEffect);
    }

    public static boolean hasCodingEffectDiff(final CodingEffect effect1, final CodingEffect effect2)
    {
        if((effect1 == UNDEFINED || effect1 == NONE) && (effect2 == UNDEFINED || effect2 == NONE))
            return false;

        return effect1 != effect2;
    }
}
