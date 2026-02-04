package com.hartwig.hmftools.sage.quality;

public class JitterVariableParams
{
    public final double Optimal4;
    public final double Optimal5;
    public final double Optimal6;
    public final double Gradient;

    public JitterVariableParams(final double optimal4, final double optimal5, final double optimal6, final double gradient)
    {
        Optimal4 = optimal4;
        Optimal5 = optimal5;
        Optimal6 = optimal6;
        Gradient = gradient;
    }
}
