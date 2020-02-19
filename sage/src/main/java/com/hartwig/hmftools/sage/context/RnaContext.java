package com.hartwig.hmftools.sage.context;

public class RnaContext {

    public static final RnaContext EMPTY = new RnaContext(0, 0, 0);

    private final int refSupport;
    private final int altSupport;
    private final int depth;

    public RnaContext(final AltContext altContext) {
        this(altContext.primaryReadContext().refSupport(),
                altContext.primaryReadContext().altSupport(),
                altContext.primaryReadContext().depth());
    }


    public RnaContext(final int refSupport, final int altSupport, final int depth) {
        this.refSupport = refSupport;
        this.altSupport = altSupport;
        this.depth = depth;
    }

    public int supportRef() {
        return refSupport;
    }

    public int supportAlt() {
        return altSupport;
    }

    public int depth() {
        return depth;
    }
}
