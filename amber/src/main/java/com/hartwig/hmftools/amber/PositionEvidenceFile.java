package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

public final class PositionEvidenceFile
{
    private static final String TUMOR_EXTENSION = ".amber.tumor.raw.tsv.gz";

    public static String generateTumorDataFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + TUMOR_EXTENSION;
    }

    private enum Columns
    {
        Chromosome,
        Position,
        Ref,
        Alt,
        ReadDepth,
        IndelCount,
        RefCount,
        AltCount,
        BaseQualFiltered,
        MapQualFiltered,
        SeqTechFiltered;
    }

    public static void write(final String filename, final List<PositionEvidence> bafs) throws IOException
    {
        BufferedWriter writer = createBufferedWriter(filename);

        StringJoiner sj = new StringJoiner(TSV_DELIM);

        Arrays.stream(Columns.values()).forEach(x -> sj.add(x.toString()));
        writer.write(sj.toString());
        writer.newLine();

        for(PositionEvidence positionEvidence : bafs)
        {
            StringJoiner rowSj = new StringJoiner(TSV_DELIM);
            rowSj.add(positionEvidence.Chromosome);
            rowSj.add(String.valueOf(positionEvidence.Position));
            rowSj.add(positionEvidence.Ref.name());
            rowSj.add(positionEvidence.Alt.name());
            rowSj.add(String.valueOf(positionEvidence.ReadDepth));
            rowSj.add(String.valueOf(positionEvidence.IndelCount));
            rowSj.add(String.valueOf(positionEvidence.RefSupport));
            rowSj.add(String.valueOf(positionEvidence.AltSupport));
            rowSj.add(String.valueOf(positionEvidence.BaseQualFiltered));
            rowSj.add(String.valueOf(positionEvidence.MapQualFiltered));
            rowSj.add(String.valueOf(positionEvidence.SeqTechFiltered));
            writer.write(rowSj.toString());
            writer.newLine();
        }

        writer.close();
    }
}
