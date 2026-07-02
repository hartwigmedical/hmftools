package com.hartwig.hmftools.orange.e2e;

import com.google.common.base.Preconditions;

public class PositionInPage
{
    final double fractionDownPage;
    final double fractionAcrossPage;

    public PositionInPage(final double fractionAcrossPage, final double fractionDownPage)
    {
        Preconditions.checkArgument(fractionAcrossPage >= 0.0);
        Preconditions.checkArgument(fractionAcrossPage <= 1.0);
        Preconditions.checkArgument(fractionDownPage >= 0.0);
        Preconditions.checkArgument(fractionDownPage <= 1.0);
        this.fractionAcrossPage = fractionAcrossPage;
        this.fractionDownPage = fractionDownPage;
    }
}
