package com.hartwig.hmftools.common.vis;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.awt.Color;

public class ColorUtil
{
    public static final Color DARK_BLUE = new Color(0, 0, 139);
    public static final Color DARK_GREEN = new Color(0, 139, 0);
    public static final Color PURPLE = new Color(128, 0, 128);

    public static Color interpolateColors(final Color colorStart, final Color colorEnd, double proportion)
    {
        proportion = max(0.0, min(1.0, proportion));

        int red = (int) Math.round((1 - proportion) * colorStart.getRed() + proportion * colorEnd.getRed());
        red = max(0, min(255, red));

        int green = (int) Math.round((1 - proportion) * colorStart.getGreen() + proportion * colorEnd.getGreen());
        green = max(0, min(255, green));

        int blue = (int) Math.round((1 - proportion) * colorStart.getBlue() + proportion * colorEnd.getBlue());
        blue = max(0, min(255, blue));

        return new Color(red, green, blue);
    }

    public static Color lighten(final Color color, double proportion)
    {
        return interpolateColors(color, Color.WHITE, proportion);
    }

    public static Color darken(final Color color, double proportion)
    {
        return interpolateColors(color, Color.BLACK, proportion);
    }
}
