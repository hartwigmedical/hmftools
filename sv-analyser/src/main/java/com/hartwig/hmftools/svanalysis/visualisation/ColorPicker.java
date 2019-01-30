package com.hartwig.hmftools.svanalysis.visualisation;

import java.awt.Color;

import org.jetbrains.annotations.NotNull;

public interface ColorPicker {

    // Dark
    Color COLOR1 = new Color(31, 120, 180);
    Color COLOR2 = new Color(227, 26, 28);
    Color COLOR3 = new Color(255, 127, 0);
    Color COLOR4 = new Color(51, 160, 44);
    Color COLOR5 = new Color(106, 61, 154);

    // Light
    Color COLOR6 = new Color(166, 206, 227);
    Color COLOR7 = new Color(251, 154, 153);
    Color COLOR8 = new Color(253, 191, 111);
    Color COLOR9 = new Color(178, 223, 138);
    Color COLOR10 = new Color(202, 178, 214);
    Color COLOR11 = new Color(255, 255, 153);

    Color[] COLOURS = new Color[] { COLOR1, COLOR3, COLOR5, COLOR6, COLOR7, COLOR8, COLOR9, COLOR10, COLOR11 };

    @NotNull
    default String color(final int clusterId, final int chainId) {
        return chainId < COLOURS.length ?  toString(COLOURS[chainId]) : "color=black";
    }

    @NotNull
    static String toString(@NotNull final  Color color) {
        return "color=(" + color.getRed() + "," + color.getGreen() + "," + color.getBlue() + ")";
    }

    @NotNull
    static String simpleSvColor(@NotNull final String type) {
        switch (type) {
            case "DEL" : return toString(COLOR2);
            case "DUP" : return toString(COLOR4);
        }

        return toString(COLOR4);
    }


}
