package com.hartwig.hmftools.common.variant.impact;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public enum VariantEffect
{
    STOP_GAINED("stop_gained", 60),
    STOP_LOST("stop_lost", 60),
    START_LOST("start_lost", 60),
    FRAMESHIFT("frameshift_variant", 60),

    SPLICE_ACCEPTOR("splice_acceptor_variant", 50),
    SPLICE_DONOR("splice_donor_variant", 50),

    INFRAME_INSERTION("inframe_insertion", 40),
    INFRAME_DELETION("inframe_deletion", 40),
    MISSENSE("missense_variant", 40),
    PHASED_INFRAME_INSERTION("phased_inframe_insertion", 40),
    PHASED_INFRAME_DELETION("phased_inframe_deletion", 40),

    SYNONYMOUS("synonymous_variant", 30),

    INTRONIC("intron_variant", 10),
    FIVE_PRIME_UTR("5_prime_UTR_variant", 10),
    THREE_PRIME_UTR("3_prime_UTR_variant", 10),
    UPSTREAM_GENE("upstream_gene_variant", 10),
    NON_CODING_TRANSCRIPT("non_coding_transcript_exon_variant", 10),

    OTHER("other", 0);

    private final String mEffect;
    private final int mRank;

    public static final String VARIANT_EFFECTS_DELIM = "&";

    VariantEffect(final String effectName, final int rank)
    {
        mEffect = effectName;
        mRank = rank;
    }

    public String effect() { return mEffect; }
    public int rank() { return mRank; }

    public boolean isEffectOf(@NotNull final String name) { return name.equals(mEffect); }

    public static List<VariantEffect> convertFromEffects(@NotNull final List<String> effects)
    {
        final List<VariantEffect> variantEffects = Lists.newArrayList();
        effects.forEach(x -> variantEffects.add(fromEffect(x)));
        return variantEffects;
    }

    public static VariantEffect fromEffect(@NotNull final String name)
    {
        for(final VariantEffect variantEffect : VariantEffect.values())
        {
            if(variantEffect.isEffectOf(name))
                return variantEffect;
        }

        return VariantEffect.OTHER;
    }

    public static String effectsToString(final List<VariantEffect> variantEffects)
    {
        return effectsToString(variantEffects, VARIANT_EFFECTS_DELIM);
    }

    public static String effectsToString(final List<VariantEffect> variantEffects, final String delim)
    {
        StringJoiner sj = new StringJoiner(delim);
        variantEffects.forEach(x -> sj.add(x.effect()));
        return sj.toString();
    }
}
