package com.hartwig.hmftools.linx.fusion;

import java.util.List;

import com.google.common.collect.Lists;

public class FusionParameters
{
    public boolean CheckInvalidReasons;
    public List<String> InvalidReasons;
    public boolean AllowExonSkipping;
    public boolean RequirePhaseMatch;

    public FusionParameters()
    {
        CheckInvalidReasons = false;
        InvalidReasons = Lists.newArrayList();
        AllowExonSkipping = true;
        RequirePhaseMatch = false;
    }

}
