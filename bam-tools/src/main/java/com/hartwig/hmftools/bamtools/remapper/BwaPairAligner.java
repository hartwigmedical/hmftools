package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.umccr.java.hellbender.utils.bwa.BwaMemAligner;
import org.umccr.java.hellbender.utils.bwa.BwaMemAlignment;
import org.umccr.java.hellbender.utils.bwa.BwaMemIndex;

public class BwaPairAligner implements PairAligner
{
    private final BwaMemAligner mAligner;

    public BwaPairAligner(BwaMemIndex index)
    {
        mAligner = new BwaMemAligner(index);
        mAligner.alignPairs();
        mAligner.dontInferPairEndStats();
        mAligner.setBandwidthOption(31);
    }

    @Override
    public ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignSequences(final byte[] bases1, final byte[] bases2)
    {
        List<List<BwaMemAlignment>> rawResults = mAligner.alignSeqs(List.of(bases1, bases2));
        return ImmutablePair.of(rawResults.get(0), rawResults.get(1));
    }
}
