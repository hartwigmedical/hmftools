package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;

import htsjdk.samtools.SAMRecord;

public interface HlaRecordAligner
{
    List<SAMRecord> alignPair(RecordPair pair);
}
