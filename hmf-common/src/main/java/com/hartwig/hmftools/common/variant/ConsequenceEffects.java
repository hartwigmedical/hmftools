package com.hartwig.hmftools.common.variant;

import java.util.List;

import com.google.common.collect.Lists;

public final class ConsequenceEffects
{
    public static final String EFFECTS_SEPARATOR = "&";

    public static final String FIVE_PRIME_UTR_EFFECT = "5_prime_UTR_variant";
    public static final String THREE_PRIME_UTR_EFFECT = "3_prime_UTR_variant";
    public static final String SPLICE_ACCEPTOR_EFFECT = "splice_acceptor_variant";
    public static final String SPLICE_DONOR_EFFECT = "splice_donor_variant";
    public static final String SPLICE_REGION_EFFECT = "splice_region_variant";
    public static final String INTRON_VARIANT_EFFECT = "intron_variant";

    public static List<String> toEffects(final String effectString)
    {
        return Lists.newArrayList(effectString.split(EFFECTS_SEPARATOR));
    }

    public static String addEffect(final String effect, final String effects, boolean atStart)
    {
        return atStart ? effect + EFFECTS_SEPARATOR + effects : effects + EFFECTS_SEPARATOR + effect;
    }
}
