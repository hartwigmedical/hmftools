package com.hartwig.hmftools.isofox.novel.cohort;

import java.util.List;

import com.google.common.collect.Lists;

public class RecurrentVariant
{
    public final String Chromsome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public final String Key;

    public final List<String> SampleIds;

    public RecurrentVariant(final String chromsome, final int position, final String ref, final String alt)
    {
        Chromsome = chromsome;
        Position = position;
        Ref = ref;
        Alt = alt;

        Key = formKey(chromsome, position, ref, alt);
        SampleIds = Lists.newArrayList();
    }

    public boolean matches(final RecurrentVariant other)
    {
        return Chromsome.equals(other.Chromsome) && Position == other.Position && Ref.equals(other.Ref) && Alt.equals(other.Alt);
    }

    public boolean matches(final String chromsome, final int position, final String ref, final String alt)
    {
        return Chromsome.equals(chromsome) && Position == position && Ref.equals(ref) && Alt.equals(alt);
    }

    public static String formKey(final String chromsome, final int position, final String ref, final String alt)
    {
        return String.format("%s-%d-%s-%s", chromsome, position, ref, alt);
    }
}
