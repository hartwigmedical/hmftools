package com.hartwig.hmftools.sage.sam;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public interface CigarHandler {

    default void handleLeftSoftClip(@NotNull final SAMRecord record, final int cigarIndex, @NotNull final CigarElement element) {
    }

    default void handleRightSoftClip(@NotNull final SAMRecord record, final int cigarIndex, @NotNull final CigarElement element,
            int readIndex, int refPosition) {
    }

    default void handleAlignment(@NotNull final SAMRecord record, final int cigarIndex, @NotNull final CigarElement element, int readIndex,
            int refPosition) {
    }

    // Note that for insert, readIndex and refPosition are BEFORE the event
    default void handleInsert(@NotNull final SAMRecord record, final int cigarIndex, @NotNull final CigarElement element, int readIndex,
            int refPosition) {
    }

    // Note that for delete, readIndex and refPosition are BEFORE the event
    default void handleDelete(@NotNull final SAMRecord record, final int cigarIndex, @NotNull final CigarElement element, int readIndex,
            int refPosition) {
    }

    // Note that for skipped, readIndex and refPosition are BEFORE the event
    default void handleSkippedReference(@NotNull final SAMRecord record, final int cigarIndex, @NotNull final CigarElement element,
            int readIndex, int refPosition) {
    }
}
