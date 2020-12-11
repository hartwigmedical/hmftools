package com.hartwig.hmftools.common.neo;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWNSTREAM;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UPSTREAM;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class NeoEpitopeFusion
{
    public final String[] GeneIds;
    public final String[] GeneNames;
    public final String[] Chromosomes;
    public final int[] Positions;
    public final byte[] Orientations;
    public final int[] SvIds;
    public final String InsertSequence;

    public static final String DELIMITER = ",";

    public NeoEpitopeFusion(
            final String geneIdUp, final String geneNameUp, final String chrUp, int posUp, byte orientUp, int svIdUp,
            final String geneIdDown, final String geneNameDown, final String chrDown, int posDown, byte orientDown, int svIdDown,
            final String insertSequence)
    {
        GeneIds = new String[] { geneIdUp, geneIdDown };
        GeneNames = new String[] { geneNameUp, geneNameDown };
        Chromosomes = new String[] { chrUp, chrDown };
        Positions = new int[] { posUp, posDown };
        Orientations = new byte[] { orientUp, orientDown };
        SvIds = new int[] { svIdUp, svIdDown };
        InsertSequence = insertSequence;
    }

    private static final String FILE_EXTENSION = ".linx.neo_epitope.tsv";
    public static final String NE_SAMPLE_ID = "sampleId";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<NeoEpitopeFusion> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<NeoEpitopeFusion> fusions) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(fusions));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<NeoEpitopeFusion> fusions)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        fusions.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<NeoEpitopeFusion> fromLines(@NotNull List<String> lines)
    {
        if(lines.isEmpty())
            return Lists.newArrayList();

        final String header = lines.get(0);

        boolean skipSampleId = header.startsWith(NE_SAMPLE_ID);
        lines.remove(0);

        return lines.stream().map(x -> NeoEpitopeFusion.fromString(x, skipSampleId)).collect(toList());
    }

    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("geneIdUp")
                .add("geneNameUp")
                .add("chromosomeUp")
                .add("positionUp")
                .add("orientationUp")
                .add("svIdUp")
                .add("geneIdDown")
                .add("geneNameDown")
                .add("chromosomeDown")
                .add("positionDown")
                .add("orientationDown")
                .add("svIdDown")
                .add("insertSeq")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final NeoEpitopeFusion fusion)
    {
        StringJoiner sj = new StringJoiner(DELIMITER);

        for(int fs = FS_UPSTREAM; fs <= FS_DOWNSTREAM; ++fs)
        {
            sj.add(String.valueOf(fusion.GeneIds[fs]));
            sj.add(String.valueOf(fusion.GeneNames[fs]));
            sj.add(String.valueOf(fusion.Chromosomes[fs]));
            sj.add(String.valueOf(fusion.Positions[fs]));
            sj.add(String.valueOf(fusion.Orientations[fs]));
            sj.add(String.valueOf(fusion.SvIds[fs]));
        }

        sj.add(fusion.InsertSequence);

        return sj.toString();
    }

    @NotNull
    public static NeoEpitopeFusion fromString(@NotNull final String data, boolean skipSampleId)
    {
        final String[] values = data.split(DELIMITER, -1);

        int index = 0;

        if(skipSampleId)
            ++index;

        return new NeoEpitopeFusion(
                values[index++], values[index++], values[index++],
                Integer.parseInt(values[index++]), Byte.parseByte(values[index++]), Integer.parseInt(values[index++]),
                values[index++], values[index++], values[index++],
                Integer.parseInt(values[index++]), Byte.parseByte(values[index++]), Integer.parseInt(values[index++]), values[index++]);
    }
}
