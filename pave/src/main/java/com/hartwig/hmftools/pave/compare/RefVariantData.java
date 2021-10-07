package com.hartwig.hmftools.pave.compare;

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
    public final String HgvsCodingImpact;
    public final String HgvsProteinImpact;

    public final CodingEffect WorstCodingEffect;
    public final int GenesAffected;

    public final String Microhomology;
    public final String RepeatSequence;
    public final boolean PhasedInframeIndel;
    public final boolean Reported;

    // repeatCount, localPhaseSet, localRealignmentSet, phasedInframeIndel, reported

    public RefVariantData(
            final String chromosome, final int position, final String ref, final String alt, final VariantType type,
            final String gene, final String canonicalEffect, final CodingEffect canonicalCodingEffect,
            final CodingEffect worstCodingEffect, final int genesAffected,
            final String hgvsCodingImpact, final String hgvsProteinImpact,
            final String microhomology, final String repeatSequence, boolean phasedInframeIndel, boolean reported)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
        Gene = gene;

        CanonicalEffect = canonicalEffect;
        CanonicalCodingEffect = canonicalCodingEffect;
        WorstCodingEffect = worstCodingEffect;
        GenesAffected = genesAffected;

        HgvsCodingImpact = hgvsCodingImpact;
        HgvsProteinImpact = hgvsProteinImpact;
        Microhomology = microhomology;
        RepeatSequence = repeatSequence;
        PhasedInframeIndel = phasedInframeIndel;
        Reported = reported;
    }

    public static RefVariantData fromSomatic(final SomaticVariant variant)
    {
        return new RefVariantData(
                variant.chromosome(), (int)variant.position(), variant.ref(), variant.alt(), variant.type(), variant.gene(),
                variant.canonicalEffect(), variant.canonicalCodingEffect(), variant.worstCodingEffect(),
                variant.genesAffected(),  variant.canonicalHgvsCodingImpact(), variant.canonicalHgvsProteinImpact(),
                variant.microhomology(), variant.repeatSequence(), variant.phasedInframeIndelIdentifier() > 0,
                variant.reported());
    }

    public String toString()
    {
        return String.format("pos(%s:%d) variant(%s: %s>%s) canon(%s: %s) worst(%s)",
                Chromosome, Position, Type, Ref, Alt, CanonicalCodingEffect, CanonicalEffect,
                WorstCodingEffect);
    }

    public static boolean hasCodingEffectDiff(final CodingEffect effect1, final CodingEffect effect2)
    {
        if((effect1 == UNDEFINED || effect1 == NONE) && (effect2 == UNDEFINED || effect2 == NONE))
            return false;

        return effect1 != effect2;
    }
}
