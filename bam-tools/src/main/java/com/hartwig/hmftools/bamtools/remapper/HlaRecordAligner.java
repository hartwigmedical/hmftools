package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public interface HlaRecordAligner
{
    @NotNull List<SAMRecord> alignPair(@NotNull RecordPair pair);
}
