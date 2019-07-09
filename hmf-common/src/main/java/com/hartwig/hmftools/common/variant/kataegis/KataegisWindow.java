package com.hartwig.hmftools.common.variant.kataegis;

import htsjdk.variant.variantcontext.VariantContext;

public class KataegisWindow {

    private final String contig;
    private final long start;

    private long end;
    private int count;

    public KataegisWindow(final VariantContext context) {
        this.contig = context.getContig();
        this.start = context.getStart();
        this.end = this.start;
        this.count = 0;
    }

    public KataegisWindow(final KataegisWindow window) {
        this.contig = window.contig;
        this.start = window.start;
        this.end = window.end;
        this.count = window.count;
    }

    public void add(final VariantContext context) {
        this.count++;
        this.end = context.getStart();
    }

    public int count() {
        return count;
    }

    public long end() {
        return end;
    }

    public boolean isViable(int minCount, long maxAverageDistance) {
        return count() >= minCount && averageDistance() <= maxAverageDistance;
    }

    public long averageDistance() {
        return count == 0 || count == 1 ? 0 : (end - start) / (count - 1);
    }

}
