package com.hartwig.hmftools.common.variant.enrich;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.utils.ArrayDeck;
import com.hartwig.hmftools.common.variant.NearPonFilteredIndel;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;

public class NearPonFilteredIndelEnrichment implements VariantContextEnrichment {

    private static final int BUFFER_SIZE = 500;
    private static final String NEAR_INDEL_PON_FILTER = "NEAR_INDEL_PON";
    private static final String NEAR_INDEL_PON_DESCRIPTION = "Indel is near a PON filtered indel";

    private final Consumer<VariantContext> consumer;
    private final ArrayDeck<VariantContext> deque = new ArrayDeck<>(BUFFER_SIZE + 1);

    public NearPonFilteredIndelEnrichment(@NotNull final Consumer<VariantContext> consumer) {
        this.consumer = consumer;
    }

    @Override
    public void accept(@NotNull final VariantContext context) {

        // We need to look at variants on both sides of victim so we need to process the variants up to midpoint (but only once)
        if (deque.size() == BUFFER_SIZE - 1) {
            for (int i = 0; i < BUFFER_SIZE / 2; i++) {
                checkAndUpdate(i);
            }
        }

        // Process variant at midpoint and flush oldest
        if (deque.size() >= BUFFER_SIZE) {
            checkAndUpdate(BUFFER_SIZE / 2);
            consumer.accept(deque.pop());
        }

        deque.add(context);
    }

    private void checkAndUpdate(int index) {
        if (NearPonFilteredIndel.isNearPonFilteredIndel(index, deque)) {
            deque.get(index).getCommonInfo().addFilter(NEAR_INDEL_PON_FILTER);
        }
    }

    @Override
    public void flush() {
        for (int i = 0; i < deque.size(); i++) {
            checkAndUpdate(i);
        }

        while (!deque.isEmpty()) {
            consumer.accept(deque.pop());
        }
    }

    @NotNull
    @Override
    public VCFHeader enrichHeader(@NotNull final VCFHeader template) {
        template.addMetaDataLine(new VCFFilterHeaderLine(NEAR_INDEL_PON_FILTER, NEAR_INDEL_PON_DESCRIPTION));
        return template;
    }
}
