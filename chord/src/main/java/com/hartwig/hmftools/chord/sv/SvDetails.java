package com.hartwig.hmftools.chord.sv;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.utils.file.FileWriterUtils;

public class SvDetails
{
    public static BufferedWriter initializeWriter(String path) throws IOException
    {
        BufferedWriter writer = FileWriterUtils.createBufferedWriter(path, false);

        String header = String.join(
                TSV_DELIM,
                "Id",

                "RefChromosome",
                "RefPosition",
                "RefOrientation",
                "RefString",

                "AltChromosome",
                "AltPosition",
                "AltOrientation",
                "AltString",

                "Type",
                "Length"
        );

        writer.write(header);
        writer.newLine();

        return writer;
    }

    public static void writeLine(BufferedWriter writer, StructuralVariant sv) throws IOException
    {
        String line = String.join(
                TSV_DELIM,
                sv.Id,

                sv.RefChromosome.toString(),
                String.valueOf(sv.RefPosition),
                sv.RefOrientation.toString(),
                sv.Context.getReference().getBaseString(),

                sv.AltChromosome.toString(),
                String.valueOf(sv.AltPosition),
                sv.AltOrientation.toString(),
                sv.Context.getAlternateAllele(0).getDisplayString(),

                sv.Type.toString(),
                String.valueOf(sv.Length)
        );

        writer.write(line);
        writer.newLine();
    }
}
