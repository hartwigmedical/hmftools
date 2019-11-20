package com.hartwig.hmftools.sage.phase;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.read.ReadContextCounter;
import com.hartwig.hmftools.sage.variant.SageVariant;
import com.hartwig.hmftools.sage.variant.SageVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;

class MnvFactory {

    private static final Logger LOGGER = LogManager.getLogger(MnvFactory.class);

    @NotNull
    private final IndexedFastaSequenceFile reference;
    private final SageVariantFactory sageVariantFactory;

    MnvFactory(@NotNull final IndexedFastaSequenceFile reference, final SageVariantFactory sageVariantFactory) {
        this.reference = reference;
        this.sageVariantFactory = sageVariantFactory;
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

        SageVariant result = sageVariantFactory.create(normal, alts);
        result.localPhaseSet(left.localPhaseSet());
        result.synthetic(true);

        return result;

    }

    @NotNull
    private AltContext merge(@NotNull final VariantHotspot variant, @NotNull final AltContext left, @NotNull final AltContext right) {
        final ReadContextCounter counter = new ReadContextCounter(variant, left.primaryReadContext(), right.primaryReadContext());
        final AltContext result = new AltContext(left.refContext(), variant.ref(), variant.alt());
        result.setPrimaryReadContext(counter);

        return result;
    }

    @NotNull
    private VariantHotspot createMnv(@NotNull final AltContext left, @NotNull final AltContext right) {
        int mnvLength = (int) (right.position() - left.position() + 1);
        int additionalLength = mnvLength - left.alt().length();

        try {
            final String alt = left.alt() + right.primaryReadContext().readContext().mnvAdditionalAlt(additionalLength);
            final String ref = reference.getSubsequenceAt(left.chromosome(), left.position(), right.position()).getBaseString();
            return ImmutableVariantHotspotImpl.builder().from(left).ref(ref).alt(alt).build();
        } catch (Exception e) {
            LOGGER.error("Unable to merge {}:{} with {}", left.chromosome(), left.position(), right.position());
            throw e;
        }

    }

}
