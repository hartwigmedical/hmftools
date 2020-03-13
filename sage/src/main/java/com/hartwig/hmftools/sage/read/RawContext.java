package com.hartwig.hmftools.sage.read;

public class RawContext {
    private final int readIndex;
    private final boolean readIndexInDelete;
    private final boolean altSupport;
    private final boolean refSupport;
    private final boolean depthSupport;

    public RawContext(final int readIndex, final boolean readIndexInDelete, final boolean altSupport, final boolean refSupport, final boolean depthSupport) {
        this.readIndex = readIndex;
        this.readIndexInDelete = readIndexInDelete;
        this.altSupport = altSupport;
        this.refSupport = refSupport;
        this.depthSupport = depthSupport;
    }

    public int readIndex() {
        return readIndex;
    }

    public boolean isReadIndexInDelete() {
        return readIndexInDelete;
    }

    public boolean isAltSupport() {
        return altSupport;
    }

    public boolean isRefSupport() {
        return refSupport;
    }

    public boolean isDepthSupport() {
        return depthSupport;
    }
}
