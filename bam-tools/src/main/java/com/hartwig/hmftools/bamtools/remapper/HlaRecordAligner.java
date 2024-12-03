package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public interface HlaRecordAligner
{
    static SAMRecord createRemappedRecord(@NotNull final SAMRecord record, @NotNull final BwaMemAlignment alignment)
    {
        SAMRecord remappedRecord = record.deepCopy();
        remappedRecord.setReferenceIndex(alignment.getRefId());
        remappedRecord.setAlignmentStart(alignment.getRefStart());
        remappedRecord.setCigarString(alignment.getCigar());
        remappedRecord.setMappingQuality(alignment.getMapQual());
//        remappedRecord.setFlags(alignment.getSamFlag());
        return remappedRecord;
    }

    @NotNull List<SAMRecord> alignPair(@NotNull RecordPair pair);
    @NotNull List<SAMRecord> alignRecord(@NotNull SAMRecord record);
}
