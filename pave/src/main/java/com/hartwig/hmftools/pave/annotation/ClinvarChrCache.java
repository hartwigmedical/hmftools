package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.pave.annotation.ClinvarAnnotation.CLNSIG;
import static com.hartwig.hmftools.pave.annotation.ClinvarAnnotation.CLNSIGCONF;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.StringCache;
import com.hartwig.hmftools.pave.VariantData;

public class ClinvarChrCache
{
    public final String Chromosome;

    private int mCurrentIndex;
    private final List<ClinvarEntry> mEntries;
    private final StringCache mStringCache;

    public ClinvarChrCache(final String chromosome, final StringCache stringCache)
    {
        Chromosome = chromosome;
        mCurrentIndex = 0;
        mEntries = Lists.newArrayList();
        mStringCache = stringCache;
    }

    public void addEntry(final int position, final String ref, final String alt, final String significance, final String conflict)
    {
        mEntries.add(new ClinvarEntry(
                position, mStringCache.intern(ref), mStringCache.intern(alt), stripBrackets(significance), stripBrackets(conflict)));
    }

    public void clear() { mEntries.clear(); }

    private static String stripBrackets(final String clinvarStr)
    {
        return clinvarStr.replaceAll(
                "\\[", "").replaceAll("\\]", "").replaceAll(" ", "");
    }

    public void annotateVariant(final VariantData variant)
    {
        if(mEntries.isEmpty() || mEntries.get(mEntries.size() - 1).Position < variant.Position)
            return;

        int position = variant.Position;
        int firstPosMatchIndex = -1;

        for(; mCurrentIndex < mEntries.size(); ++mCurrentIndex)
        {
            ClinvarEntry entry = mEntries.get(mCurrentIndex);

            if(entry.Position < position)
                continue;

            if(entry.Position > position)
                break;

            if(firstPosMatchIndex == -1)
                firstPosMatchIndex = mCurrentIndex;

            if(entry.matches(variant))
            {
                variant.context().getCommonInfo().putAttribute(CLNSIG, entry.Significance);

                if(!entry.Conflict.isEmpty())
                    variant.context().getCommonInfo().putAttribute(CLNSIGCONF, entry.Conflict);

                break;
            }
        }

        // move the index back to the prior position or the first at this position
        if(firstPosMatchIndex >= 0)
            mCurrentIndex = firstPosMatchIndex;
        else if(mCurrentIndex > 0)
            --mCurrentIndex;
    }

    private class ClinvarEntry
    {
        public final int Position;
        public final String Ref;
        public final String Alt;
        public final String Significance;
        public final String Conflict;

        public ClinvarEntry(final int position, final String ref, final String alt, final String significance, final String conflict)
        {
            Position = position;
            Ref = ref;
            Alt = alt;
            Significance = significance;
            Conflict = conflict;
        }

        public boolean matches(final VariantData variant)
        {
            return variant.Position == Position && variant.Ref.equals(Ref) && variant.Alt.equals(Alt);
        }

        public String toString() { return String.format("%d %s>%s details(%s - %s)",
                Position, Ref, Alt, Significance, Conflict); }
    }
}
