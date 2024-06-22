package com.hartwig.hmftools.redux;

import static com.hartwig.hmftools.redux.merge.BamMerger.formSequenceIntervals;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.redux.merge.SequenceInfo;

import org.junit.Test;

import htsjdk.samtools.SAMSequenceRecord;

public class BamMergeTest
{
    @Test
    public void testSequenceInfos()
    {

        List<SAMSequenceRecord> sequences = Lists.newArrayList();

        sequences.add(new SAMSequenceRecord("1", 10000));
        sequences.add(new SAMSequenceRecord("2", 10000));
        sequences.add(new SAMSequenceRecord("3", 10000));

        for(int i = 0; i < sequences.size(); ++i)
        {
            sequences.get(i).setSequenceIndex(i);
        }

        List<SequenceInfo> sequenceInfoList = formSequenceIntervals(sequences, "bam_file", 4);
        assertEquals(4, sequenceInfoList.size());

        SequenceInfo seqInfo = sequenceInfoList.get(0);
        assertTrue(hasInterval(seqInfo, 0, 1, 7500));

        seqInfo = sequenceInfoList.get(1);
        assertTrue(hasInterval(seqInfo, 0, 7501, 10000));
        assertTrue(hasInterval(seqInfo, 1, 1, 5000));

        seqInfo = sequenceInfoList.get(2);
        assertTrue(hasInterval(seqInfo, 1, 5001, 10000));
        assertTrue(hasInterval(seqInfo, 2, 1, 2500));

        seqInfo = sequenceInfoList.get(3);
        assertTrue(hasInterval(seqInfo, 2, 2501, 10000));

        sequenceInfoList = formSequenceIntervals(sequences, "bam_file", 1);
        assertEquals(1, sequenceInfoList.size());

        seqInfo = sequenceInfoList.get(0);
        assertTrue(hasInterval(seqInfo, 0, 1, 10000));
        assertTrue(hasInterval(seqInfo, 1, 1, 10000));
        assertTrue(hasInterval(seqInfo, 2, 1, 10000));

        sequenceInfoList = formSequenceIntervals(sequences, "bam_file", 2);
        assertEquals(2, sequenceInfoList.size());

        seqInfo = sequenceInfoList.get(0);
        assertTrue(hasInterval(seqInfo, 0, 1, 10000));
        assertTrue(hasInterval(seqInfo, 1, 1, 5000));

        seqInfo = sequenceInfoList.get(1);
        assertTrue(hasInterval(seqInfo, 1, 5001, 10000));
        assertTrue(hasInterval(seqInfo, 2, 1, 10000));
    }

    private static boolean hasInterval(final SequenceInfo seqInfo, int seqIndex, int regionStart, int regionEnd)
    {
        return seqInfo.Intervals.stream().anyMatch(x -> x.referenceIndex == seqIndex && x.start == regionStart && x.end == regionEnd);
    }
}
