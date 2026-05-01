package com.hartwig.hmftools.purple.plot;

import static com.hartwig.hmftools.purple.plot.CircosFileWriter.circosContig;
import static com.hartwig.hmftools.purple.plot.CircosFileWriter.writeCircosFile;

import java.io.IOException;
import java.util.Collection;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.variant.VariantContextDecorator;

public final class CircosIndelWriter
{
    public static void writePositions(final String filePath, Collection<VariantContextDecorator> values)
            throws IOException
    {
        writeCircosFile(filePath, values, CircosIndelWriter::transformPosition);
    }

    private static String transformPosition(final VariantContextDecorator variant)
    {
        return new StringJoiner("\t").add(circosContig(variant.chromosome()))
                .add(String.valueOf(variant.position()))
                .add(String.valueOf(variant.position()))
                .add("1") // VAF has no impact on the line (since type=heatmap)
                .add("color=" + color(variant))
                .toString();
    }

    private static String color(final VariantContextDecorator variant)
    {
        // Inserts: GREEN
        if(variant.alt().length() > variant.ref().length())
        {
            return "green";
        }

        // DEL with repeats >=4: ORANGE
        if(variant.repeatCount() >= 4)
        {
            // del repeat >= 4
            return "orange";
        }

        // DEL with microhomology > 2: RED
        if(variant.microhomology().length() > 2)
        {
            return "vdred";
        }


        // other DELs: YELLOW
        return "vdyellow";
    }
}
