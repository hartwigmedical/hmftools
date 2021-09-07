package com.hartwig.hmftools.linx.ext_compare;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

public class CfLink
{
    public final CfSvData[] SVs;
    public final int[] LinkIndex;

    public CfLink(final CfSvData var1, final CfSvData var2, final int linkIndex1, final int linkIndex2)
    {
        SVs = new CfSvData[] { var1, var2 };
        LinkIndex = new int[] { linkIndex1, linkIndex2 };
    }

    public int length() { return abs(SVs[SE_START].Positions[LinkIndex[SE_START]] - SVs[SE_END].Positions[LinkIndex[SE_END]]); }

    public CfSvData getOtherSv(final CfSvData var) { return SVs[SE_START] == var ? SVs[SE_END] : SVs[SE_START]; }

}
