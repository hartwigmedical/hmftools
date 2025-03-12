package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;

import htsjdk.samtools.SAMRecord;

public interface SingleRecordAligner
{
    List<SAMRecord> alignSequence(SAMRecord original);
}
