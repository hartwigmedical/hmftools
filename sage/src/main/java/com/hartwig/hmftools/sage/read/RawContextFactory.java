package com.hartwig.hmftools.sage.read;

import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.sam.CigarTraversal;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class RawContextFactory {

    private final VariantHotspot variant;

    public RawContextFactory(final VariantHotspot variant) {
        this.variant = variant;
    }

    @NotNull
    public RawContext create(@NotNull final SAMRecord record) {
        RawContextCigarHandler handler = new RawContextCigarHandler(variant);
        CigarTraversal.traverseCigar(record, handler);
        RawContext result = handler.result();
        return result == null ? new RawContext(-1, false, false, false, false) : result;

    }

}
