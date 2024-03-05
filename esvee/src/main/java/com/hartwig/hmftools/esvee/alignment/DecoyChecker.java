package com.hartwig.hmftools.esvee.alignment;

import static com.hartwig.hmftools.esvee.SvConstants.DECOY_MAX_MISMATCHES;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

public class DecoyChecker
{
    private final BwaMemAligner mAligner;
    private int mSequenceCount;
    private int mSequenceMatched;

    public DecoyChecker(final String decoyGenome)
    {
        mSequenceCount = 0;
        mSequenceMatched = 0;

        if(decoyGenome != null && Files.exists(Paths.get(decoyGenome)))
        {
            BwaMemIndex index = new BwaMemIndex(decoyGenome);
            mAligner = new BwaMemAligner(index);
        }
        else
        {
            mAligner = null;
        }
    }

    public boolean enabled() { return mAligner != null;}

    public int sequenceCount() { return mSequenceCount; }
    public int sequenceMatched() { return mSequenceMatched; };

    public boolean matchesDecoy(final String sequence)
    {
        ++mSequenceCount;

        List<BwaMemAlignment> alignmentResults = mAligner.alignSeqs(List.of(sequence.getBytes())).get(0);

        for(BwaMemAlignment alignment : alignmentResults)
        {
            if(alignment.getNMismatches() <= DECOY_MAX_MISMATCHES)
            {
                ++mSequenceMatched;
                return true;
            }
        }

        return false;
    }
}
