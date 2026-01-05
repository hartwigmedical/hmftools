package com.hartwig.hmftools.sage.seqtech;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.redux.BaseQualAdjustment.isLowBaseQual;

import com.hartwig.hmftools.sage.common.VariantReadContext;

import htsjdk.samtools.SAMRecord;

public class IlluminaUtils
{
    public static boolean lowQualInVariantCore(int variantReadIndex, final VariantReadContext readContext, final SAMRecord read)
    {
        if(readContext.variant().isIndel())
        {
            if(readContext.MaxRepeat != null && readContext.MaxRepeat.Count >= 15)
                return false;

            int coreIndexStart = max(variantReadIndex - readContext.leftCoreLength(), 0);
            int coreIndexEnd = min(variantReadIndex + readContext.rightCoreLength(), read.getBaseQualities().length - 1);

            for(int i = coreIndexStart; i <= coreIndexEnd; ++i)
            {
                if(isLowBaseQual(read.getBaseQualities()[i]))
                    return true;
            }
        }
        else
        {
            for(int i = variantReadIndex; i <= variantReadIndex + readContext.variant().altLength() - 1; ++i)
            {
                if(isLowBaseQual(read.getBaseQualities()[i]))
                    return true;
            }
        }

        return false;
    }
}
