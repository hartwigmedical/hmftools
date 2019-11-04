package com.hartwig.hmftools.sage.count;

import java.util.ArrayDeque;
import java.util.function.Consumer;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.SAMRecord;

@Deprecated
public class SAMRecordQueue implements Consumer<SAMRecord> {

    private final Consumer<SAMRecord> consumer;
    private final ArrayDeque<SAMRecord> samRecords = new ArrayDeque<>();

    public SAMRecordQueue(final Consumer<SAMRecord> consumer) {
        this.consumer = consumer;
    }

    public void flush() {
        SAMRecord last = samRecords.peekLast();
        if (last != null) {
            flushUntil(last.getAlignmentEnd() + 1);
        }
    }

    public void flushUntil(long alignmentEnd) {
        while (true) {
            @Nullable
            SAMRecord oldest = samRecords.peekFirst();
            if (oldest != null && oldest.getAlignmentEnd() <= alignmentEnd) {
                consumer.accept(oldest);
                samRecords.removeLast();
            } else {
                return;
            }
        }
    }

    @Override
    public void accept(final SAMRecord record) {
        samRecords.addLast(record);
    }
}
