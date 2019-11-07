package com.hartwig.hmftools.sage.phase;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class PhasingMnv {

    @NotNull
    final IndexedFastaSequenceFile reference;

    public PhasingMnv(@NotNull final IndexedFastaSequenceFile reference) {
        this.reference = reference;
    }

    @NotNull
    public SageVariant createMNV(@NotNull final SageVariant left, @NotNull final SageVariant right) {

        final VariantHotspot variant = createMnv(left.normal(), right.normal());
        final AltContext normal = merge(variant, left.normal(), right.normal());

        final List<AltContext> alts = Lists.newArrayList();
        for (int i = 0; i < left.tumorAltContexts().size(); i++) {
            final AltContext leftTumor = left.tumorAltContexts().get(i);
            final AltContext rightTumor = right.tumorAltContexts().get(i);
            alts.add(merge(variant, leftTumor, rightTumor));
        }

        return new SageVariant(normal, alts);

    }

    @NotNull
    private AltContext merge(@NotNull final VariantHotspot variant, @NotNull final AltContext left, @NotNull final AltContext right) {

        final ReadContextCounter counter = new ReadContextCounter(variant, left.primaryReadContext(), right.primaryReadContext());

        final AltContext result = new AltContext(left.sample(), variant);
        result.setPrimaryReadContext(counter);
        return result;
    }

    @NotNull
    private VariantHotspot createMnv(@NotNull final AltContext left, @NotNull final AltContext right) {
        int length = (int) (right.position() - left.position() + 1);

        final String alt = left.primaryReadContext().readContext().alt(length);
        final String ref = reference.getSubsequenceAt(left.chromosome(), left.position(), right.position()).getBaseString();

        return ImmutableVariantHotspotImpl.builder().from(left).ref(ref).alt(alt).build();
    }

}
