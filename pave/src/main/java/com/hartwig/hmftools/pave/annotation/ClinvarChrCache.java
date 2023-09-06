package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.pave.annotation.ClinvarAnnotation.CLNSIG;
import static com.hartwig.hmftools.pave.annotation.ClinvarAnnotation.CLNSIGCONF;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.pave.VariantData;

public class ClinvarChrCache
{
    public final String Chromosome;

    private int mCurrentIndex;
    private final List<ClinvarEntry> mEntries;

    public ClinvarChrCache(final String chromosome)
    {
        Chromosome = chromosome;
        mCurrentIndex = 0;
        mEntries = Lists.newArrayList();
    }

    public void addEntry(final int position, final String ref, final String alt, final String significance, final String conflict)
    {
        mEntries.add(new ClinvarEntry(position, ref, alt, stripBrackets(significance), stripBrackets(conflict)));
    }

    public void clear() { mEntries.clear(); }

    private static String stripBrackets(final String clinvarStr)
    {
        return clinvarStr.replaceAll(
                "\\[", "").replaceAll("\\]", "").replaceAll(" ", "");
    }

    public void annotateVariant(final VariantData variant)
    {
        // Clinvar entries for both v37 and v38 do not have the chromosome prefix

        if(mEntries.isEmpty() || mEntries.get(mEntries.size() - 1).Position < variant.Position)
            return;

        int position = variant.Position;

        for(; mCurrentIndex < mEntries.size(); ++mCurrentIndex)
        {
            ClinvarEntry entry = mEntries.get(mCurrentIndex);

            if(entry.Position < position)
                continue;

            if(entry.Position > position)
            {
                if(mCurrentIndex > 0)
                    --mCurrentIndex;

                return;
            }

            if(entry.matches(variant))
            {
                variant.context().getCommonInfo().putAttribute(CLNSIG, entry.Significance);

                if(!entry.Conflict.isEmpty())
                    variant.context().getCommonInfo().putAttribute(CLNSIGCONF, entry.Conflict);
            }

            return;
        }
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
