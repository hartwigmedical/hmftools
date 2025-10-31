package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.esvee.common.SvConstants.isUltima;

public final class SeqTechUtils
{
    public static final int READ_MISMATCH_MEDIUM_REPEAT_COUNT_ULTIMA = 7;
    public static final int READ_MISMATCH_LONG_REPEAT_COUNT_ULTIMA = 11;

    public static void setSeqTechSpecifics()
    {
        if(isUltima())
        {
            AssemblyConstants.READ_MISMATCH_MEDIUM_REPEAT_COUNT = READ_MISMATCH_MEDIUM_REPEAT_COUNT_ULTIMA;
            AssemblyConstants.READ_MISMATCH_LONG_REPEAT_COUNT = READ_MISMATCH_LONG_REPEAT_COUNT_ULTIMA;
        }
    }
}
