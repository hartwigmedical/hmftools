package com.hartwig.hmftools.sage.dedup;

import static com.hartwig.hmftools.sage.dedup.DedupIndel.dedupIndels;
import static com.hartwig.hmftools.sage.dedup.DedupMatching.dedupMatchingVariants;
import static com.hartwig.hmftools.sage.dedup.DedupSnvMnv.dedupMnvOverlaps;
import static com.hartwig.hmftools.sage.dedup.DedupSnvMnv.dedupMnvSnvs;

import java.util.List;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.common.SageVariant;

public class VariantDeduper
{
    private final DedupMixedGermlineSomatic mDedupMixedGermlineSomatic;

    public VariantDeduper(final List<TranscriptData> transcripts)
    {
        mDedupMixedGermlineSomatic = new DedupMixedGermlineSomatic(transcripts);
    }

    public void processVariants(final List<SageVariant> variants)
    {
        dedupMnvOverlaps(variants);

        mDedupMixedGermlineSomatic.dedupVariants(variants);

        dedupMnvSnvs(variants);

        dedupIndels(variants);

        dedupMatchingVariants(variants);
    }

    private void acceptDedupedVariant(final SageVariant variant) {}

    public static boolean longerContainsShorter(final SageVariant shorter, final SageVariant longer)
    {
        return longerContainsShorter(shorter.variant(), longer.variant());
    }

    public static boolean longerContainsShorter(final VariantHotspot shorter, final VariantHotspot longer)
    {
        int longerStart = longer.position();
        int longerEnd = longer.end();

        int shorterStart = shorter.position();
        int shorterEnd = shorter.end();

        if(shorterStart < longerStart || shorterEnd > longerEnd)
            return false;

        final String shorterAlt = shorter.alt();

        int offset = shorterStart - longerStart;
        final String longerAlt = new String(longer.alt().getBytes(), offset, shorter.alt().length());
        return shorterAlt.equals(longerAlt);
    }

}
