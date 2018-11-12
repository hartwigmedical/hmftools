package com.hartwig.hmftools.common.hotspot;

import java.util.function.Consumer;

import com.hartwig.hmftools.common.pileup.Pileup;

class MNVEvidence implements Consumer<Pileup> {

    private int reads = Integer.MAX_VALUE;
    private int refCount = Integer.MAX_VALUE;
    private int altCount = Integer.MAX_VALUE;
    private int score = Integer.MAX_VALUE;

    private final VariantHotspot mnv;

    MNVEvidence(final VariantHotspot mvn) {
        this.mnv = mvn;
    }

    public int reads() {
        return reads;
    }

    public int altCount() {
        return altCount == Integer.MAX_VALUE ? 0 : altCount;
    }

    public int refCount() {
        return refCount == Integer.MAX_VALUE ? 0 : refCount;
    }

    public int score() {
        return score;
    }

    @Override
    public void accept(final Pileup pileup) {
        int positionOffset = (int) (pileup.position() - mnv.position());

        char alt = mnv.alt().charAt(positionOffset);
        refCount = Math.min(refCount, pileup.referenceCount());
        altCount = Math.min(altCount, pileup.mismatchCount(alt));
        score = Math.min(score, pileup.mismatchScore(alt));
        reads = Math.min(reads, pileup.readCount());
    }
}
