package com.hartwig.hmftools.common.samtools;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public interface CigarHandler {

    default void handleLeftSoftClip(final SAMRecord record, final CigarElement element) {}

    default void handleRightSoftClip(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    default void handleAlignment(final SAMRecord record, final CigarElement element, boolean beforeIndel, int readIndex, int refPosition) {}

    // Note that for insert, readIndex and refPosition are BEFORE the event
    default void handleInsert(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    // Note that for delete, readIndex and refPosition are BEFORE the event
    default void handleDelete(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}

    // Note that for skipped, readIndex and refPosition are BEFORE the event
    default void handleSkippedReference(final SAMRecord record, final CigarElement element, int readIndex, int refPosition) {}
}
