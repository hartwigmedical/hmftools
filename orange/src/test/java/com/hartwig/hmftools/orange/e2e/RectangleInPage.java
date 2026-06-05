package com.hartwig.hmftools.orange.e2e;

import com.google.common.base.Preconditions;

import org.junit.Assert;

public class RectangleInPage
{
    final PositionInPage topLeft;
    final PositionInPage bottomRight;

    public RectangleInPage(double topLeftX, double topLeftY, double bottomRightX, double bottomRightY)
    {
        this(new PositionInPage(topLeftX, topLeftY), new PositionInPage(bottomRightX, bottomRightY));
    }

    public RectangleInPage(final PositionInPage topLeft, final PositionInPage bottomRight)
    {
        Preconditions.checkArgument(topLeft.fractionAcrossPage < bottomRight.fractionAcrossPage);
        Preconditions.checkArgument(topLeft.fractionDownPage < bottomRight.fractionDownPage);
        this.bottomRight = bottomRight;
        this.topLeft = topLeft;
    }

    public boolean contains(final PositionInPage position)
    {
        return position.fractionAcrossPage >= topLeft.fractionAcrossPage &&
                position.fractionAcrossPage <= bottomRight.fractionAcrossPage &&
                position.fractionDownPage >= topLeft.fractionDownPage &&
                position.fractionDownPage <= bottomRight.fractionDownPage;
    }

    public boolean contains(RectangleInPage other)
    {
        return contains(other.topLeft) && contains(other.bottomRight);
    }
}
