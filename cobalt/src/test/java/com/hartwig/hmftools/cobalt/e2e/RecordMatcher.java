package com.hartwig.hmftools.cobalt.e2e;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import htsjdk.samtools.SAMRecord;

public class RecordMatcher
{
    private final Random mRandom = new Random();
    private final double FractionOfRecordsToKeep;
    private final Map<String, SAMRecord> mRecordsByName = new HashMap<>();

    private long Processed = 0;
    private long Retained = 0;

    public RecordMatcher(final double fractionOfRecordsToKeep)
    {
        FractionOfRecordsToKeep = fractionOfRecordsToKeep;
    }

    public List<SAMRecord> process(final SAMRecord record)
    {
        ++Processed;

        if(mRecordsByName.containsKey(record.getReadName()))
        {
            SAMRecord match = mRecordsByName.remove(record.getReadName());
            if (mRandom.nextDouble() <= FractionOfRecordsToKeep)
            {
                Retained += 2;
                return List.of(record, match);
            }
            else
            {
                return List.of();
            }
        }
        else
        {
            mRecordsByName.put(record.getReadName(), record);
            return List.of();
        }
    }

    public long processed()
    {
        return Processed;
    }

    public long retained()
    {
        return Retained;
    }
}
