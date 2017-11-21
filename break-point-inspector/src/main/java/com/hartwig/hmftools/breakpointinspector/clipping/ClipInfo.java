package com.hartwig.hmftools.breakpointinspector.clipping;

import com.hartwig.hmftools.breakpointinspector.Location;

import htsjdk.samtools.SAMRecord;

public class ClipInfo {

    final SAMRecord Record;

    ClipInfo(final SAMRecord record) {
        this.Record = record;
    }

    Location Alignment;
    int Length = 0;
    String Sequence = "";
    boolean HardClipped = false;
    boolean Left = false;
    boolean Right = false;
}
