package com.hartwig.hmftools.purple.plot;

import static com.hartwig.hmftools.purple.plot.CircosFileWriter.circosContig;
import static com.hartwig.hmftools.purple.plot.CircosFileWriter.writeCircosFile;

import java.io.IOException;
import java.util.Collection;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.variant.VariantContextDecorator;

public final class CircosINDELWriter
{
    public static void writePositions(final String filePath, Collection<VariantContextDecorator> values)
            throws IOException
    {
        writeCircosFile(filePath, values, CircosINDELWriter::transformPosition);
    }

    private static String transformPosition(final VariantContextDecorator position)
    {
        return new StringJoiner("\t").add(circosContig(position.chromosome()))
                .add(String.valueOf(position.position()))
                .add(String.valueOf(position.position()))
                .add(String.valueOf(position.adjustedVaf()))
                .add("fill_color=" + color(position))
                .toString();
    }

    //    DEL REPEAT>=4  ORANGE
    //    DEL MH>2 RED
    //    OTHER DEL YELLOW
    //    INS  GREEN
    private static String color(final VariantContextDecorator variant)
    {
        if(variant.alt().length() > variant.ref().length())
        {
            // insertion
            return "green";
        }

        if(variant.repeatCount() >= 4)
        {
            // del repeat >= 4
            return "orange";
        }

        if(!variant.microhomology().isEmpty())
        {
            // microhomology
            return "vdred";
        }

        return "vdyellow";
    }
}
