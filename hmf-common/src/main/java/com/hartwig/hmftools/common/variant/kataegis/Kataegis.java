package com.hartwig.hmftools.common.variant.kataegis;

import java.util.function.Consumer;
import java.util.function.Predicate;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class Kataegis implements Consumer<VariantContext> {

    private final Predicate<VariantContext> forwardPredicate;
    private final Predicate<VariantContext> reversePredicate;

    private final KataegisQueue forwardDetector;
    private final KataegisQueue reverseDetector;

    public Kataegis(final Consumer<VariantContext> consumer) {

        reverseDetector = new KataegisQueue(KataegisStatus.FWD, x -> true, consumer);
        forwardDetector = new KataegisQueue(KataegisStatus.FWD, x -> true, reverseDetector);

        forwardPredicate = x -> true;
        reversePredicate = x -> true;
    }

    @Override
    public void accept(@NotNull final VariantContext context) {
        forwardDetector.accept(context);
    }

    public void flush() {
        forwardDetector.flush();
        reverseDetector.flush();
    }
}
