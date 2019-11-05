package com.hartwig.hmftools.sage.sam;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public interface CigarHandler {

    void handleAlignment(@NotNull final SAMRecord record, @NotNull final CigarElement element, int readIndex, int refPosition);

    // Note that for insert, readIndex and refPosition are BEFORE the event
    void handleInsert(@NotNull final SAMRecord record, @NotNull final CigarElement element, int readIndex, int refPosition);

    // Note that for delete, readIndex and refPosition are BEFORE the event
    void handleDelete(@NotNull final SAMRecord record, @NotNull final CigarElement element, int readIndex, int refPosition);

}
