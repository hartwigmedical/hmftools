package com.hartwig.hmftools.qsee.prep.category;

import static org.junit.Assert.assertEquals;

import static com.hartwig.hmftools.qsee.prep.category.ReduxMsIndelErrorPrep.DELETIONS_KEY;
import static com.hartwig.hmftools.qsee.prep.category.ReduxMsIndelErrorPrep.INSERTIONS_KEY;

import java.util.List;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.bam.ConsensusType;
import com.hartwig.hmftools.common.redux.JitterTableRow;

import org.junit.Test;

import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.prep.category.ReduxMsIndelErrorPrep;

public class ReduxMsIndelErrorPrepTest
{
    @Test
    public void canAggregateByJitter()
    {
        JitterTableRow input = createJitterTableRow(
                new int[] { -3, -2, -1, 0, 1, 2, 3 },
                new int[] {  1,  2,  5, 0, 6, 3, 1 }
        );

        JitterTableRow output = ReduxMsIndelErrorPrep.aggregateByJitter(List.of(input)).get(0);

        int actualJitterCountDeletions = output.jitterCounts().get(DELETIONS_KEY);
        int actualJitterCountInsertions = output.jitterCounts().get(INSERTIONS_KEY);

        assertEquals(8, actualJitterCountDeletions);
        assertEquals(10, actualJitterCountInsertions);
    }

    @Test
    public void canAggregateByRepeatType()
    {
        List<JitterTableRow> table = List.of(
                createJitterTableRow(0, "A/T", ConsensusType.SINGLE, 1, 1, 2),
                createJitterTableRow(0, "A/T", ConsensusType.SINGLE, 1, 1, 2),

                // Change ref units
                createJitterTableRow(1, "A/T", ConsensusType.SINGLE, 2, 2, 4),
                createJitterTableRow(1, "A/T", ConsensusType.SINGLE, 2, 2, 4),

                // Change consensus type
                // 1 nucleotide repeat
                createJitterTableRow(1, "A/T", ConsensusType.NONE, 3, 3, 6),
                createJitterTableRow(1, "A/T", ConsensusType.NONE, 3, 3, 6),

                createJitterTableRow(1, "C/G", ConsensusType.NONE, 4, 4, 8),
                createJitterTableRow(1, "C/G", ConsensusType.NONE, 4, 4, 8),

                // 2 nucleotide repeat
                createJitterTableRow(1, "NN/NN", ConsensusType.NONE, 5, 5, 10),
                createJitterTableRow(1, "MM/MM/MM/MM", ConsensusType.NONE, 5, 5, 10),

                // >=3 nucleotide repeat
                createJitterTableRow(1, "1bp repeat", ConsensusType.NONE, 6, 6, 12),
                createJitterTableRow(1, "2bp repeat", ConsensusType.NONE, 6, 6, 12),
                createJitterTableRow(1, "3bp repeat", ConsensusType.NONE, 6, 6, 12)
        );

        List<JitterTableRow> expectedOutput = List.of(
                createJitterTableRow(0, "A/T repeat", ConsensusType.SINGLE, 2, 2, 4),
                createJitterTableRow(1, "A/T repeat", ConsensusType.SINGLE, 4, 4, 8),
                createJitterTableRow(1, "A/T repeat", ConsensusType.NONE, 6, 6, 12),
                createJitterTableRow(1, "C/G repeat", ConsensusType.NONE, 8, 8, 16),
                createJitterTableRow(1, "2bp repeat", ConsensusType.NONE, 10, 10, 20),
                createJitterTableRow(1, "≥3bp repeat", ConsensusType.NONE, 18, 18, 36)
        );

        List<JitterTableRow> actualOutput = ReduxMsIndelErrorPrep.aggregateByRepeatType(table);

        assertEquals(expectedOutput.size(), actualOutput.size());

        for(int i = 0; i < expectedOutput.size(); ++i)
        {
            JitterTableRow expectedRow = expectedOutput.get(i);
            JitterTableRow actualRow = actualOutput.get(i);

            assertEquals(expectedRow.refNumUnits(), actualRow.refNumUnits());
            assertEquals(expectedRow.getRepeatUnit(), actualRow.getRepeatUnit());
            assertEquals(expectedRow.getConsensusType(), actualRow.getConsensusType());
            assertEquals(expectedRow.totalReadCount(), actualRow.totalReadCount());
            assertEquals(expectedRow.jitterCounts(), actualRow.jitterCounts());
        }
    }

    @Test
    public void canCalculatePhredScores()
    {
        List<JitterTableRow> table = List.of(
                createJitterTableRow(4, "A/T repeat", ConsensusType.NONE, 6, 4, 100)
        );

        List<Feature> actualOutput = ReduxMsIndelErrorPrep.calcPhredScores(table);

        Feature microsatelliteIndelErrorRate = actualOutput.get(0);
        Feature microsatelliteIndelErrorBias = actualOutput.get(1);

        assertEquals(FeatureType.MS_INDEL_ERROR_RATES, microsatelliteIndelErrorRate.key().type());
        assertEquals(10.0, microsatelliteIndelErrorRate.value(), 0.01);

        assertEquals(FeatureType.MS_INDEL_ERROR_BIAS, microsatelliteIndelErrorBias.key().type());
        assertEquals(-1.76, microsatelliteIndelErrorBias.value(), 0.01);
    }

    private static JitterTableRow createJitterTableRow(int[] jitters, int[] counts)
    {
        JitterTableRow row = new JitterTableRow(0, "", null);

        row.setTotalReadCount(IntStream.of(counts).sum());

        for(int i = 0; i < jitters.length; ++i)
        {
            row.setJitterReadCount(jitters[i], counts[i]);
        }

        return row;
    }

    private static JitterTableRow createJitterTableRow(int refNumUnits, String repeatUnit, ConsensusType consensusType,
            int deletionReadCount, int insertionCount, int totalReadCount)
    {
        JitterTableRow row = new JitterTableRow(refNumUnits, repeatUnit, consensusType);

        row.setJitterReadCount(DELETIONS_KEY, deletionReadCount);
        row.setJitterReadCount(INSERTIONS_KEY, insertionCount);
        row.setTotalReadCount(totalReadCount);

        return row;
    }
}
