package com.hartwig.hmftools.bamtools.remapper;

import static com.hartwig.hmftools.bamtools.remapper.RemapperTestBase.bwa;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.codon.Nucleotides;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.umccr.java.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;

import htsjdk.samtools.SAMRecord;

/**
 * This class implements Aligner by hard-coding the alignments for the hla alt records
 * in the test data file tiny.sam to the values returned by a BwaAligner instance that uses
 * the hg38_no_alts reference sequence.
 */
class HlaAligner implements PairAligner
{

    private final List<SAMRecord> mRecords;
    private final Map<Integer, List<BwaMemAlignment>> mRecordIndexToAlignments = new HashMap<>();

    HlaAligner(List<SAMRecord> records)
    {
        mRecords = records;

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:1516:12156:4225", records.get(4).getReadName()); // Sanity check
        // 81,0,194358862,194358978,35,151,60,0,116,66,35S116M,116,null,5,29945439,0
        // 161,5,29945439,29945590,0,151,60,0,151,116,151M,151,null,0,194358862,0
        mRecordIndexToAlignments.put(4, List.of(bwa("81,0,194358862,194358978,35,151,60,0,116,66,35S116M,116,null,5,29945439,0")));

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:1516:12156:4225", records.get(7).getReadName()); // Sanity check
        mRecordIndexToAlignments.put(7, List.of(bwa("161,5,29945439,29945590,0,151,60,0,151,116,151M,151,null,0,194358862,0")));

        // 81,5,31354760,31354911,0,151,60,2,145,101,151M,21G128T0,null,5,31354347,-564
        // 161,5,31354347,31354419,0,72,60,2,62,19,72M79S,4C45C21,null,5,31354760,564
        // 2225,1,32916241,32916273,42,74,13,0,32,28,42S32M77S,32,null,5,31354760,0
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:2559:3224:21292", records.get(8).getReadName()); // Sanity check
        Assert.assertTrue(records.get(8).getReadNegativeStrandFlag()); // Sanity check
        mRecordIndexToAlignments.put(8, List.of(
                bwa("81,5,31354760,31354911,0,151,60,2,145,101,151M,21G128T0,null,5,31354347,-564")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:1:2559:3224:21292", records.get(9).getReadName()); // Sanity check
        Assert.assertFalse(records.get(9).getReadNegativeStrandFlag()); // Sanity check
        mRecordIndexToAlignments.put(9, List.of(
                bwa("161,5,31354347,31354419,0,72,60,2,62,19,72M79S,4C45C21,null,5,31354760,564"),
                bwa("2225,1,32916241,32916273,42,74,13,0,32,28,42S32M77S,32,null,5,31354760,0")
        ));

        // 97,5,31355729,31355880,0,151,60,0,151,19,151M,151,null,5,31356297,719
        // 145,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,5,31355729,-719
        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1446:18213:29684", records.get(10).getReadName()); // Sanity check
        Assert.assertFalse(records.get(10).getReadNegativeStrandFlag()); // Sanity check
        mRecordIndexToAlignments.put(10, List.of(
                bwa("97,5,31355729,31355880,0,151,60,0,151,19,151M,151,null,5,31356297,719")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:1:1446:18213:29684", records.get(11).getReadName()); // Sanity check
        mRecordIndexToAlignments.put(11, List.of(
                bwa("145,5,31356297,31356448,0,151,60,1,146,108,151M,76C74,null,5,31355729,-719")
        ));

        // 97,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,5,31354513,289
        // 2161,1,32916486,32916522,67,103,0,0,36,36,67S36M48S,36,null,5,31354513,0
        // 145,5,31354513,31354664,0,151,60,1,146,84,151M,111C39,null,5,31354375,-289
        Assert.assertEquals("A00624:8:HHKYHDSXX:4:1304:7075:14998", records.get(12).getReadName()); // Sanity check
        mRecordIndexToAlignments.put(12, List.of(
                bwa("145,5,31354513,31354664,0,151,60,1,146,84,151M,111C39,null,5,31354375,-289")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:4:1304:7075:14998", records.get(13).getReadName()); // Sanity check
        mRecordIndexToAlignments.put(13, List.of(
                bwa("97,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,5,31354513,289"),
                bwa("2161,1,32916486,32916522,67,103,0,0,36,36,67S36M48S,36,null,5,31354513,0")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:2276:27299:7623", records.get(14).getReadName()); // Sanity check
        mRecordIndexToAlignments.put(14, List.of(
                bwa("97,5,31356204,31356355,0,151,19,4,131,121,151M,21T0C19C79G28", "chr6,+31271111,151M,6;", "5,31271513,-84542")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:2:2276:27299:7623", records.get(15).getReadName()); // Sanity check
        mRecordIndexToAlignments.put(15, List.of(
                bwa("145,5,31271513,31271664,0,151,18,6,121,111,151M,5G8C0T28G42T51C11", "chr6,-31356606,48M3I100M,10;", "5,31356204,84542")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:4:2543:5737:19288", records.get(16).getReadName()); // Sanity check
        mRecordIndexToAlignments.put(16, List.of(
                bwa("121,5,29943927,29944078,0,151,60,0,151,73,151M,151,null,-1,-1,0")
        ));

        Assert.assertEquals("A00624:8:HHKYHDSXX:4:2543:5737:19288", records.get(17).getReadName()); // Sanity check
        mRecordIndexToAlignments.put(17, List.of(
                bwa("181,-1,-1,-1,-1,-1,0,0,0,0,\"\",null,null,5,29943927,0")
        ));
    }

    public List<BwaMemAlignment> alignSequence(byte[] bases)
    {
        for(int index = 0; index < mRecords.size(); index++)
        {
            final byte[] recordBases = mRecords.get(index).getReadBases();
            if(Arrays.equals(bases, recordBases))
            {
                return mRecordIndexToAlignments.get(index);
            }
            if(Arrays.equals(bases, Nucleotides.reverseComplementBases(recordBases)))
            {
                return mRecordIndexToAlignments.get(index);
            }
        }
        throw new IllegalArgumentException("No alignment found for " + new String(bases));
    }

    @Override
    public ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignSequences(final byte[] bases1, final byte[] bases2)
    {
        List<BwaMemAlignment> resultLeft = alignSequence(bases1);
        List<BwaMemAlignment> resultRight = alignSequence(bases2);
        return new ImmutablePair<>(resultLeft, resultRight);
    }
}
