package com.hartwig.hmftools.purple.plot;

import static com.hartwig.hmftools.purple.plot.CircosFileWriter.circosContig;
import static com.hartwig.hmftools.purple.plot.CircosFileWriter.writeCircosFile;

import java.io.IOException;
import java.util.Collection;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;

public final class CircosSnvWriter
{
    private static final String PINK = "233,187,184";
    private static final String BLUE = "20,176,239";
    private static final String BLACK = "6,8,9";
    private static final String RED = "224,7,20";
    private static final String GREY = "191,190,191";
    private static final String GREEN = "144,202,75";

    public static void writePositions(final String filePath, Collection<VariantContextDecorator> values) throws IOException
    {
        writeCircosFile(filePath, values, CircosSnvWriter::transformPosition);
    }

    private static String transformPosition(final VariantContextDecorator variant)
    {
        return String.join("\t", circosContig(variant.chromosome()),
                String.valueOf(variant.position()),
                String.valueOf(variant.position()),
                String.valueOf(variant.adjustedVaf()),
                String.join(",", "color=" + color(variant), "glyph=" + glyph(), "glyph_size=" + glyphSize()));
    }

    private static String color(final VariantContextDecorator variant)
    {
        if(signature("C", "A", variant))
        {
            return BLUE;
        }
        if(signature("C", "G", variant))
        {
            return BLACK;
        }
        if(signature("C", "T", variant))
        {
            return RED;
        }
        if(signature("T", "A", variant))
        {
            return GREY;
        }
        if(signature("T", "C", variant))
        {
            return GREEN;
        }
        if(signature("T", "G", variant))
        {
            return PINK;
        }

        return "purple";
    }

    private static boolean signature(final String ref, final String alt, final VariantContextDecorator variant)
    {
        return (variant.ref().equals(ref) && variant.alt().equals(alt)) || (variant.ref().equals(inverse(ref)) && variant.alt()
                .equals(inverse(alt)));
    }

    private static String inverse(final String base)
    {
        return Nucleotides.swapDnaBase(base);
    }

    private static String glyph()
    {
        return "circle";
    }

    private static int glyphSize()
    {
        return 10;
    }
}