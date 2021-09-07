package com.hartwig.hmftools.sage.read;

import com.hartwig.hmftools.common.samtools.CigarTraversal;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class RawContextFactory
{
    private final VariantHotspot mVariant;

    private static final RawContext DUMMY = RawContext.inSoftClip(-1);

    public RawContextFactory(final VariantHotspot variant)
    {
        mVariant = variant;
    }

    @NotNull
    public RawContext create(final int maxSkippedReferenceRegions, @NotNull final SAMRecord record)
    {
        RawContextCigarHandler handler = new RawContextCigarHandler(maxSkippedReferenceRegions, mVariant);
        CigarTraversal.traverseCigar(record, handler);
        RawContext result = handler.result();
        return result == null ? DUMMY : result;
    }
}
