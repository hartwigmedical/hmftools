package com.hartwig.hmftools.breakpointinspector.clipping;

import com.hartwig.hmftools.breakpointinspector.Location;

import htsjdk.samtools.SAMRecord;

public class ClipInfo {

    final SAMRecord record;

    ClipInfo(final SAMRecord record) {
        this.record = record;
    }

    Location alignment;
    int length = 0;
    String sequence = "";
    boolean hardClipped = false;
    boolean left = false;
    boolean right = false;
}
