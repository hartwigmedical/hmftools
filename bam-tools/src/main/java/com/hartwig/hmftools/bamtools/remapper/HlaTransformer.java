package com.hartwig.hmftools.bamtools.remapper;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.bam.SamRecordUtils;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class HlaTransformer
{
    static boolean hasSomeHlaReference(SAMRecord record)
    {
        return hasAltReference(record) || mateHasAltReference(record);
    }

    static boolean hasAltReference(SAMRecord record)
    {
        return isHlaAltReference(record.getReferenceName());
    }

    static boolean mateHasAltReference(SAMRecord record)
    {
        return isHlaAltReference(record.getMateReferenceName());
    }

    static boolean isHlaAltReference(String referenceName)
    {
        return referenceName.toLowerCase().startsWith("hla-");
    }

    private @NotNull
    final HlaRecordAligner aligner;
    private final Map<String, SAMRecord> recordsByName = new HashMap<>();
    private int NumberProcessed = 0;

    public HlaTransformer(@NotNull final HlaRecordAligner aligner)
    {
        this.aligner = aligner;
    }

    public @NotNull List<SAMRecord> process(final @NotNull SAMRecord record)
    {
        if(!hasSomeHlaReference(record))
        {
            return List.of(record);
        }

        // Ignore supplementary records.
        if (record.isSecondaryOrSupplementary())
        {
            return List.of();
        }
        NumberProcessed++;
        if (recordsByName.containsKey(record.getReadName())) {
            SAMRecord match = recordsByName.remove(record.getReadName());
            RecordPair pair = pair(match, record);
            return aligner.alignPair(pair);
        } else {
            recordsByName.put(record.getReadName(), record);
            return List.of();
        }
    }

    public @NotNull List<SAMRecord> unmatchedRecords()
    {
        return new ArrayList<>(recordsByName.values());
    }

    public int numberOfHlaRecordsProcessed()
    {
        return NumberProcessed;
    }

    private RecordPair pair(final SAMRecord s, final SAMRecord r)
    {
        if (SamRecordUtils.firstInPair(s))
        {
            return new RecordPair(s, r);
        }
        return new RecordPair(r, s);
    }
}
