package com.hartwig.hmftools.sage.pipeline;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.context.AltContext;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.ReferenceSequenceFile;

public interface MnvPipeline {

    Logger LOGGER = LogManager.getLogger(MnvPipeline.class);

    @NotNull
    VariantHotspot combined(@NotNull final SageVariant left, @NotNull final SageVariant right);

    @Nullable
    SageVariant mnv(int lps, @NotNull final VariantHotspot mnv);


    @NotNull
    default  VariantHotspot combined(@NotNull final ReferenceSequenceFile refGenome, @NotNull final AltContext left, @NotNull final AltContext right) {
        int mnvLength = (int) (right.position() - left.position() + 1);
        int additionalLength = mnvLength - left.alt().length();

        try {
            final String alt = left.alt() + right.primaryReadContext().readContext().mnvAdditionalAlt(additionalLength);
            final String ref = refGenome.getSubsequenceAt(left.chromosome(), left.position(), right.position()).getBaseString();
            return ImmutableVariantHotspotImpl.builder().from(left).ref(ref).alt(alt).build();
        } catch (Exception e) {
            LOGGER.error("Unable to merge {}:{} with {}", left.chromosome(), left.position(), right.position());
            throw e;
        }
    }

}
