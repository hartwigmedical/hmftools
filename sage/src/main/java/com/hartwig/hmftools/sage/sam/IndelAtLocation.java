package com.hartwig.hmftools.sage.sam;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class IndelAtLocation implements CigarHandler {

    public static boolean indelAtPosition(int position, @NotNull final SAMRecord record) {
        final IndelAtLocation handler = new IndelAtLocation(position);
        CigarTraversal.traverseCigar(record, handler);
        return handler.indel();
    }

    private IndelAtLocation(final int position) {
        this.position = position;
    }

    private final int position;
    private boolean indel;

    @Override
    public void handleAlignment(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
            final int refPosition) {
    }

    @Override
    public void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
            final int refPosition) {
        indel |= refPosition == position;
    }

    @Override
    public void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement element, final int readIndex,
            final int refPosition) {
        indel |= refPosition == position;
    }

    private boolean indel() {
        return indel;
    }
}
