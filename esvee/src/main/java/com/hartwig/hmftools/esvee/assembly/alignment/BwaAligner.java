package com.hartwig.hmftools.esvee.assembly.alignment;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.BWA_PENALTY_ADJUST;
import static com.hartwig.hmftools.esvee.common.SvConstants.MIN_INDEL_LENGTH;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

public class BwaAligner implements Aligner
{
    private static final int SCORING_MATRIX_SIZE = 5;

    private final BwaMemAligner mAligner;

    public BwaAligner(final String refGenomeImageFile)
    {
        if(!refGenomeImageFile.isEmpty() && Files.exists(Paths.get(refGenomeImageFile)))
        {
            BwaMemIndex index = null;

            try
            {
                index = new BwaMemIndex(refGenomeImageFile);
            }
            catch(Exception e)
            {
                SV_LOGGER.error("failed to initialise BWA aligner: {}", e.toString());
            }

            if(index != null)
            {
                mAligner = new BwaMemAligner(index);
                mAligner.setBandwidthOption(MIN_INDEL_LENGTH - 1);

                int mismatchPenalty = mAligner.getMismatchPenaltyOption();
                mAligner.setMismatchPenaltyOption(mismatchPenalty + BWA_PENALTY_ADJUST);

                int gapOpenPenaltyDel = mAligner.getDGapOpenPenaltyOption();
                mAligner.setDGapOpenPenaltyOption(gapOpenPenaltyDel + BWA_PENALTY_ADJUST);

                int gapOpenPenaltyInsert = mAligner.getIGapOpenPenaltyOption();
                mAligner.setIGapOpenPenaltyOption(gapOpenPenaltyInsert + BWA_PENALTY_ADJUST);
                updateScoringMatrix();
            }
            else
            {
                mAligner = null;
            }
        }
        else
        {
            mAligner = null;
        }
    }

    private void updateScoringMatrix()
    {
        int matchScore = mAligner.getMatchScoreOption();
        int mismatchPenalty = mAligner.getMismatchPenaltyOption();
        byte[] scoringMatrix = new byte[SCORING_MATRIX_SIZE * SCORING_MATRIX_SIZE];
        int k = 0;

        for(int i = 0; i < SCORING_MATRIX_SIZE - 1; i++)
        {
            for(int j = 0; j < SCORING_MATRIX_SIZE - 1; j++)
                scoringMatrix[k++] = (byte) (i == j ? matchScore : -mismatchPenalty);

            scoringMatrix[k++] = -1;
        }

        for(int j = 0; j < SCORING_MATRIX_SIZE; j++)
            scoringMatrix[k++] = -1;

        mAligner.setScoringMatrixOption(scoringMatrix);
    }

    @Override
    public List<BwaMemAlignment> alignSequence(final byte[] bases)
    {
        if(mAligner == null)
            return Collections.emptyList();

        List<BwaMemAlignment> alignmentSet = mAligner.alignSeqs(List.of(bases)).get(0);

        return alignmentSet;
    }
}
