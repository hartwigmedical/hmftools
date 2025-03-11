package com.hartwig.hmftools.bamtools.remapper;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAligner;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.bwa.BwaMemIndex;

public class BwaPairAligner implements PairAligner
{
    private final BwaMemAligner mAligner;

    public BwaPairAligner(final String refGenomeImageFile)
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
                BT_LOGGER.error("failed to initialise BWA aligner: {}", e.toString());
            }

            if(index != null)
            {
                mAligner = new BwaMemAligner(index);
                mAligner.alignPairs();
                mAligner.dontInferPairEndStats();
                mAligner.setBandwidthOption(31);
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

    @Override
    public ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignSequences(final byte[] bases1, final byte[] bases2)
    {
        List<List<BwaMemAlignment>> rawResults = mAligner.alignSeqs(List.of(bases1, bases2));
        return ImmutablePair.of(rawResults.get(0), rawResults.get(1));
    }
}
