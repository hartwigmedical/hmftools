package com.hartwig.hmftools.amber;

import static java.lang.String.format;
import static java.lang.String.valueOf;

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

        if(bafs.size() > 0 && bafs.get(0).LohData != null)
        {
            sj.add("ValidBafHigh").add("ValidBafMed").add("ValidBafLow");
            sj.add("ValidRateHigh").add("ValidRateMed").add("ValidRateLow");
            sj.add("LohCandidate").add("LohRecovered");
        }

        writer.write(sj.toString());
        writer.newLine();

        for(PositionEvidence positionEvidence : bafs)
        {
            StringJoiner rowSj = new StringJoiner(TSV_DELIM);
            rowSj.add(positionEvidence.Chromosome);
            rowSj.add(valueOf(positionEvidence.Position));
            rowSj.add(positionEvidence.Ref.name());
            rowSj.add(positionEvidence.Alt.name());
            rowSj.add(valueOf(positionEvidence.ReadDepth));
            rowSj.add(valueOf(positionEvidence.IndelCount));
            rowSj.add(valueOf(positionEvidence.RefSupport));
            rowSj.add(valueOf(positionEvidence.AltSupport));
            rowSj.add(valueOf(positionEvidence.BaseQualFiltered));
            rowSj.add(valueOf(positionEvidence.MapQualFiltered));
            rowSj.add(valueOf(positionEvidence.SeqTechFiltered));

            if(positionEvidence.LohData != null)
            {
                BafLohData lohData = positionEvidence.LohData;
                rowSj.add(valueOf(lohData.ValidBafHigh));
                rowSj.add(valueOf(lohData.ValidBafMedium));
                rowSj.add(valueOf(lohData.ValidBafLow));
                rowSj.add(format("%.2f", lohData.ValidRateHigh));
                rowSj.add(format("%.2f", lohData.ValidRateMedium));
                rowSj.add(format("%.2f", lohData.ValidRateLow));
                rowSj.add(valueOf(lohData.LohCandidate));
                rowSj.add(valueOf(lohData.LohRecovered));
            }

            writer.write(rowSj.toString());
            writer.newLine();
        }

        writer.close();
    }
}
