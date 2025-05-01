package com.hartwig.hmftools.pave.annotation;

import static com.hartwig.hmftools.pave.annotation.ClinvarAnnotation.CLNSIG;
import static com.hartwig.hmftools.pave.annotation.ClinvarAnnotation.CLNSIGCONF;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.pathogenic.Pathogenicity;
import com.hartwig.hmftools.common.utils.StringCache;
import com.hartwig.hmftools.common.variant.SimpleVariant;
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

    public void addEntry(final ClinvarEntry entry) { mEntries.add(entry); }

    public void clear() { mEntries.clear(); }
    public void resetSearch() { mCurrentIndex = 0; }
    public List<ClinvarEntry> entries() { return mEntries; }
    public int entryCount() { return mEntries.size(); }

    private static String stripBrackets(final String clinvarStr)
    {
        return clinvarStr.replaceAll(
                "\\[", "").replaceAll("\\]", "").replaceAll(" ", "");
    }

    public void annotateVariant(final VariantData variant)
    {
        ClinvarEntry matchedEntry = findClinvarEntry(variant.Position, variant.Ref, variant.Alt);

        if(matchedEntry != null)
        {
            variant.context().getCommonInfo().putAttribute(CLNSIG, matchedEntry.Significance);

            if(!matchedEntry.Conflict.isEmpty())
                variant.context().getCommonInfo().putAttribute(CLNSIGCONF, matchedEntry.Conflict);
        }
    }

    private ClinvarEntry findClinvarEntry(final int position, final String ref, final String alt)
    {
        if(mEntries.isEmpty() || mEntries.get(mEntries.size() - 1).Position < position)
            return null;

        int firstPosMatchIndex = -1;

        ClinvarEntry matchedEntry = null;

        for(; mCurrentIndex < mEntries.size(); ++mCurrentIndex)
        {
            ClinvarEntry entry = mEntries.get(mCurrentIndex);

            if(entry.Position < position)
                continue;

            if(entry.Position > position)
                break;

            if(firstPosMatchIndex == -1)
                firstPosMatchIndex = mCurrentIndex;

            if(entry.matches(position, ref, alt))
            {
                matchedEntry = entry;
                break;
            }
        }

        // move the index back to the prior position or the first at this position
        if(firstPosMatchIndex >= 0)
            mCurrentIndex = firstPosMatchIndex;
        else if(mCurrentIndex > 0)
            --mCurrentIndex;

        return matchedEntry;
    }

    public Pathogenicity findPathogenicity(final SimpleVariant variant)
    {
        ClinvarEntry matchedEntry = findClinvarEntry(variant.Position, variant.Ref, variant.Alt);

        if(matchedEntry != null)
            return Pathogenicity.fromClinvarAnnotation(matchedEntry.Significance, matchedEntry.Conflict);

        return null;
    }

}
