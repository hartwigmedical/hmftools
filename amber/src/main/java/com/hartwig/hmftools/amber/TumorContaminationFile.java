package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.BaseDepthData;
import com.hartwig.hmftools.common.amber.ImmutableBaseDepthData;

import org.jetbrains.annotations.NotNull;

public final class TumorContaminationFile
{
    private static final String HEADER_PREFIX = "chr";

    private static final String AMBER_EXTENSION = ".amber.contamination.tsv";

    public static String generateContaminationFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + AMBER_EXTENSION;
    }

    public static void write(@NotNull final String filename, @NotNull final List<TumorContamination> contamination) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(contamination));
    }

    @NotNull
    static List<TumorContamination> fromLines(@NotNull List<String> lines)
    {
        final List<TumorContamination> result = Lists.newArrayList();
        for(String line : lines)
        {
            if(!line.startsWith(HEADER_PREFIX))
            {
                result.add(fromString(line));
            }
        }

        return result;
    }

    @NotNull
    private static TumorContamination fromString(@NotNull final String line)
    {
        String[] values = line.split(TSV_DELIM);

        PositionEvidence template = new PositionEvidence(values[0], Integer.parseInt(values[1]), values[2], values[3]);

        BaseDepthData normalDepth = ImmutableBaseDepthData.builder()
                .ref(BaseDepthData.Base.valueOf(template.ref()))
                .alt(BaseDepthData.Base.valueOf(template.alt()))
                .readDepth(Integer.parseInt(values[4]))
                .refSupport(Integer.parseInt(values[5]))
                .altSupport(Integer.parseInt(values[6]))
                .build();

        BaseDepthData tumorDepth = ImmutableBaseDepthData.builder()
                .ref(BaseDepthData.Base.valueOf(template.ref()))
                .alt(BaseDepthData.Base.valueOf(template.alt()))
                .readDepth(Integer.parseInt(values[7]))
                .refSupport(Integer.parseInt(values[8]))
                .altSupport(Integer.parseInt(values[9]))
                .build();

        return new TumorContamination(template.Chromosome, template.Position, normalDepth, tumorDepth);
    }

    @NotNull
    static List<String> toLines(@NotNull final List<TumorContamination> contamination)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        contamination.stream().map(TumorContaminationFile::toString).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static String header()
    {
        return new StringJoiner(TSV_DELIM, "", "")
                .add("chromosome")
                .add("position")
                .add("ref")
                .add("alt")
                .add("normalDepth")
                .add("normalRefSupport")
                .add("normalAltSupport")
                .add("tumorDepth")
                .add("tumorRefSupport")
                .add("tumorAltSupport")
                .toString();
    }

    private static String toString(final TumorContamination ratio)
    {
        return new StringJoiner(TSV_DELIM).add(ratio.chromosome())
                .add(String.valueOf(ratio.position()))
                .add(String.valueOf(ratio.Tumor.ref()))
                .add(String.valueOf(ratio.Tumor.alt()))
                //                .add(String.valueOf(ratio.Normal.readDepth()))
                .add("0")
                //                .add(String.valueOf(ratio.Normal.refSupport()))
                .add("0")
                //                .add(String.valueOf(ratio.Normal.altSupport()))
                .add("0")
                .add(String.valueOf(ratio.Tumor.readDepth()))
                .add(String.valueOf(ratio.Tumor.refSupport()))
                .add(String.valueOf(ratio.Tumor.altSupport()))
                .toString();
    }
}
