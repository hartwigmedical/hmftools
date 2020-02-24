package com.hartwig.hmftools.sage.sam;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class SkippedAtLocation implements CigarHandler {

    public static boolean inSkipReference(int maxSkipAllowed, int position, @NotNull final SAMRecord record) {
        final SkippedAtLocation handler = new SkippedAtLocation(maxSkipAllowed, position);
        CigarTraversal.traverseCigar(record, handler);
        return handler.skipped();
    }

    private SkippedAtLocation(final int maxSkip, final int position) {
        this.position = position;
        this.maxSkip = maxSkip;
    }

    private final int maxSkip;
    private final int position;
    private boolean skipped;


    @Override
    public void handleSkippedReference(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex, final int refPosition) {
        if (element.getLength() > maxSkip) {
            skipped |= (position > refPosition && position <= refPosition + element.getLength());
        }
    }

    private boolean skipped() {
        return skipped;
    }
}
