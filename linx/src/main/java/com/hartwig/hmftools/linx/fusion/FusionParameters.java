package com.hartwig.hmftools.linx.fusion;

import java.util.List;

import com.google.common.collect.Lists;

public class FusionParameters
{
    public boolean LogInvalidReasons;
    public List<String> InvalidReasons;
    public boolean AllowExonSkipping;
    public boolean RequirePhaseMatch;
    public boolean RequireUpstreamBiotypes;

    public FusionParameters()
    {
        LogInvalidReasons = false;
        InvalidReasons = Lists.newArrayList();
        AllowExonSkipping = true;
        RequirePhaseMatch = false;
        RequireUpstreamBiotypes = false;
    }

}
