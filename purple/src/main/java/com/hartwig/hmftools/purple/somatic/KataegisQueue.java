package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.KATAEGIS_FLAG;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.function.Consumer;
import java.util.function.Predicate;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class KataegisQueue {

    static final long MAX_ABS_DISTANCE = 2000;
    static final long MAX_AVG_DISTANCE = 1000;
    static final int MIN_COUNT = 3;

    private final Predicate<VariantContext> candidate;
    private final Consumer<VariantContext> consumer;
    private final Deque<VariantContext> buffer;
    private final String idPrefix;

    private int identifier = 0;

    public KataegisQueue(final String idPrefix, final Predicate<VariantContext> candidate, final Consumer<VariantContext> consumer) {
        this.candidate = candidate;
        this.consumer = consumer;
        this.idPrefix = idPrefix;
        this.buffer = new ArrayDeque<>();
    }

    public void accept(@NotNull final VariantContext context) {
        if (!buffer.isEmpty()) {
            final VariantContext previous = buffer.peekLast();
            if (context.getStart() - previous.getStart() > MAX_ABS_DISTANCE || !previous.getContig().equals(context.getContig())) {
                flush();
            }
        }

        buffer.add(context);
    }

    public void flush() {
        while (!buffer.isEmpty()) {
            processFirstContext();
        }
    }

    private void processFirstContext() {
        if (!buffer.isEmpty()) {
            final VariantContext first = buffer.peekFirst();
            if (!candidate.test(first)) {
                consumer.accept(buffer.pollFirst());
            } else {
                final KataegisWindow window = longestViableWindow(first);
                final boolean isWindowViable = window.isViable(MIN_COUNT, MAX_AVG_DISTANCE);
                if (isWindowViable) {
                    identifier++;
                }

                while (!buffer.isEmpty()) {
                    final VariantContext peek = buffer.peekFirst();
                    if (peek.getStart() > window.end()) {
                        return;
                    }

                    if (isWindowViable && candidate.test(peek)) {
                        peek.getCommonInfo().putAttribute(KATAEGIS_FLAG, idPrefix + "_" + identifier, true);
                    }

                    consumer.accept(buffer.pollFirst());
                }
            }
        }
    }

    @NotNull
    private KataegisWindow longestViableWindow(@NotNull final VariantContext first) {
        KataegisWindow result = new KataegisWindow(first);
        final KataegisWindow window = new KataegisWindow(first);

        for (VariantContext context : buffer) {
            if (context.getStart() - window.end() > MAX_ABS_DISTANCE) {
                return result;
            }

            if (candidate.test(context)) {
                window.add(context);
            }

            if (window.isViable(MIN_COUNT, MAX_AVG_DISTANCE)) {
                result = new KataegisWindow(window);
            }
        }

        return result;
    }
}
