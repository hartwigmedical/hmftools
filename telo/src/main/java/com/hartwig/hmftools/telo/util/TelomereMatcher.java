package com.hartwig.hmftools.telo.util;

// helper class to match a sequence to telomere
// and give a similarity rate

import java.util.List;

import com.hartwig.hmftools.common.aligner.SequenceAligner;
import com.hartwig.hmftools.telo.TeloConstants;
import com.hartwig.hmftools.telo.TeloUtils;

import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class TelomereMatcher
{
    private static final Logger LOGGER = LogManager.getLogger(TelomereMatcher.class);

    public static double calcCTelomereMatch(String seq)
    {
        return calcGTelomereMatch(TeloUtils.reverseComplementSequence(seq));
    }

    // calculate how much this sequence resembles telomere
    public static double calcGTelomereMatch(String seq)
    {
        // we want to allow some insertion / deletion plus a full hexamer
        int telomereTemplateLength = (int)(seq.length() * 1.2) + 6;

        // we provide a template that is 1.5x as long as the original sequence
        List<SequenceAligner.AlignOp> alignOps = SequenceAligner.alignSubsequence(seq, generateTelomereTemplate(telomereTemplateLength));

        // remove all the ops at tail that are deletes
        for (int i = alignOps.size() - 1; i >= 0; --i)
        {
            if (alignOps.get(i) == SequenceAligner.AlignOp.DELETION)
            {
                alignOps.remove(i);
            }
            else
            {
                break;
            }
        }

        // count up all the ones in the matching portion
        boolean alignStart = false;
        int numMatch = 0;
        int alignLength = 0;
        for (SequenceAligner.AlignOp alignOp : alignOps)
        {
            if (!alignStart && alignOp != SequenceAligner.AlignOp.DELETION)
            {
                alignStart = true;
            }

            if (alignStart)
            {
                ++alignLength;
                if (alignOp == SequenceAligner.AlignOp.MATCH)
                {
                    ++numMatch;
                }
            }
        }

        double matchRatio = (double)numMatch / alignLength;

        LOGGER.trace("seq({}) alignLength({}) numMatch({}) ratio({})", seq, alignLength, numMatch, matchRatio);

        return matchRatio;
    }

    private static String generateTelomereTemplate(int length)
    {
        StringBuilder b = new StringBuilder();
        int numFullHexamers = length / TeloConstants.CANONICAL_TELOMERE_SEQ.length();
        b.append(StringUtils.repeat(TeloConstants.CANONICAL_TELOMERE_SEQ, numFullHexamers));
        int residualLength = length % TeloConstants.CANONICAL_TELOMERE_SEQ.length();
        if (residualLength > 0)
        {
            b.append(TeloConstants.CANONICAL_TELOMERE_SEQ, 0, residualLength);
        }
        return b.toString();
    }
}
