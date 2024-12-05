package com.hartwig.hmftools.bamtools.remapper;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public interface HlaRecordAligner
{

    static int mergeFlags(final int original, final int fromAlignment)
    {
        Set<SAMFlag> alignmentFlags = SAMFlag.getFlags(original);
        // Remove flags that should be overridden by their value in the
        // alignment. We don't remove the supplementary flag as we never
        // re-align supplementary reads.
        alignmentFlags.remove(SAMFlag.READ_UNMAPPED);
        alignmentFlags.remove(SAMFlag.READ_REVERSE_STRAND);
        alignmentFlags.addAll(SAMFlag.getFlags(fromAlignment));
        return alignmentFlags.stream().map(SAMFlag::intValue).reduce(Integer::sum).orElse(0);
    }

    static SAMRecord createRemappedRecord(@NotNull final SAMRecord record, @NotNull final BwaMemAlignment alignment)
    {
        SAMRecord remappedRecord = record.deepCopy();
        remappedRecord.setReferenceIndex(alignment.getRefId());
        remappedRecord.setAlignmentStart(alignment.getRefStart());
        remappedRecord.setCigarString(alignment.getCigar());
        remappedRecord.setMappingQuality(alignment.getMapQual());
        remappedRecord.setFlags(mergeFlags(record.getFlags(), alignment.getSamFlag()));
        return remappedRecord;
    }

    @NotNull List<SAMRecord> alignPair(@NotNull RecordPair pair);
    @NotNull List<SAMRecord> alignRecord(@NotNull SAMRecord record);
}
