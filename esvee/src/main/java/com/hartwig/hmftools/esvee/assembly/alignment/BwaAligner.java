package com.hartwig.hmftools.esvee.assembly.alignment;

import static com.hartwig.hmftools.esvee.assembly.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.assembly.AssemblyConstants.BWA_MISMATCH_PENALTY;
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
                // mAligner.setMismatchPenaltyOption(BWA_MISMATCH_PENALTY);
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
    public List<BwaMemAlignment> alignSequence(final byte[] bases)
    {
        if(mAligner == null)
            return Collections.emptyList();

        List<BwaMemAlignment> alignmentSet = mAligner.alignSeqs(List.of(bases)).get(0);

        return alignmentSet;
    }
}
