package com.hartwig.hmftools.amber.contamination;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.common.amber.AmberBase;
import com.hartwig.hmftools.common.amber.BaseDepthData;

import org.jetbrains.annotations.NotNull;

public final class TumorContaminationFile
{
    private static final String AMBER_EXTENSION = ".amber.contamination.tsv";

    public static String generateContaminationFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + AMBER_EXTENSION;
    }

    public static void write(final String filename, final List<TumorContamination> contamination) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(contamination));
    }

    public static List<TumorContamination> read(final String filename)
    {
        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));
            return fromLines(lines);
        }
        catch(IOException e)
        {
            AMB_LOGGER.error("failed to read contamination file({}): {}", filename, e.toString());
            return null;
        }
    }

    public static List<TumorContamination> fromLines(final List<String> lines)
    {
        List<TumorContamination> contaminationEntries = Lists.newArrayList();

        if(lines.get(0).startsWith(Columns.chromosome.toString()))
            lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            String chromosome = values[Columns.chromosome.ordinal()];
            int position = Integer.parseInt(values[Columns.position.ordinal()]);
            String ref = values[Columns.ref.ordinal()];
            String alt = values[Columns.alt.ordinal()];

            BaseDepthData normalDepth = new BaseDepthData(
                    AmberBase.valueOf(ref), AmberBase.valueOf(alt),
                    Integer.parseInt(values[Columns.normalDepth.ordinal()]),
                    0,
                    Integer.parseInt(values[Columns.normalRefSupport.ordinal()]),
                    Integer.parseInt(values[Columns.normalAltSupport.ordinal()]));

            BaseDepthData tumorDepth = new BaseDepthData(
                    AmberBase.valueOf(ref), AmberBase.valueOf(alt),
                    Integer.parseInt(values[Columns.tumorDepth.ordinal()]),
                    0,
                    Integer.parseInt(values[Columns.tumorRefSupport.ordinal()]),
                    Integer.parseInt(values[Columns.tumorAltSupport.ordinal()]));

            TumorContamination tumorContamination = new TumorContamination(chromosome, position, normalDepth, tumorDepth);

            contaminationEntries.add(tumorContamination);
        }

        return contaminationEntries;
    }

    private enum Columns
    {
        chromosome,
        position,
        ref,
        alt,
        normalDepth,
        normalRefSupport,
        normalAltSupport,
        tumorDepth,
        tumorRefSupport,
        tumorAltSupport;
    }

    public static List<String> toLines(final List<TumorContamination> contamination)
    {
        List<String> lines = Lists.newArrayList();

        StringJoiner sj = new StringJoiner(TSV_DELIM);
        Arrays.stream(Columns.values()).forEach(x -> sj.add(x.toString()));
        lines.add(sj.toString());
        contamination.stream().map(TumorContaminationFile::toString).forEach(lines::add);
        return lines;
    }

    private static String toString(final TumorContamination ratio)
    {
        return new StringJoiner(TSV_DELIM).add(ratio.chromosome())
                .add(String.valueOf(ratio.position()))
                .add(String.valueOf(ratio.Tumor.ref()))
                .add(String.valueOf(ratio.Tumor.alt()))
                .add(String.valueOf(ratio.Normal.readDepth()))
                .add(String.valueOf(ratio.Normal.refSupport()))
                .add(String.valueOf(ratio.Normal.altSupport()))
                .add(String.valueOf(ratio.Tumor.readDepth()))
                .add(String.valueOf(ratio.Tumor.refSupport()))
                .add(String.valueOf(ratio.Tumor.altSupport()))
                .toString();
    }
}
