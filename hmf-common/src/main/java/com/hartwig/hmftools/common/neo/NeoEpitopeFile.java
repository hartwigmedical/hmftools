package com.hartwig.hmftools.common.neo;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFusion.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public class NeoEpitopeFile
{
    public final NeoEpitopeType VariantType;
    public final String VariantInfo;
    public final double CopyNumber;
    public final String[] GeneIds;
    public final String[] GeneNames;
    public final int[] NmdBases;
    public final String UpstreamAA;
    public final String DownstreamAA;
    public final String NovelAA;
    public final String WildtypeAA;
    public final int[] CodingBasePositions;
    public final String[] CodingBases;
    public final String[] Transcripts;

    public NeoEpitopeFile(
            final NeoEpitopeType varType, final String varInfo, final double copyNumber,
            final String geneIdUp, final String geneIdDown, final String geneNameUp, final String geneNameDown,
            final String upAA, final String downAAs, final String novelAAs, final String wildtypeAAs, int nmdBasesMin, int nmdBasesMax,
            int codingBasePosUp, int codingBasePosDown, final String upCodingBases, final String downCodingBases,
            final String transcriptsUp, final String transcriptsDown)
    {
        VariantType = varType;
        VariantInfo = varInfo;
        CopyNumber = copyNumber;
        GeneIds = new String[] { geneIdUp, geneIdDown };
        GeneNames = new String[] { geneNameUp, geneNameDown };
        NmdBases = new int[] { nmdBasesMin, nmdBasesMax };
        UpstreamAA = upAA;
        DownstreamAA = downAAs;
        NovelAA = novelAAs;
        WildtypeAA = wildtypeAAs;
        CodingBases = new String[] { upCodingBases, downCodingBases };
        CodingBasePositions = new int[] { codingBasePosUp, codingBasePosDown };
        Transcripts = new String[] { transcriptsUp, transcriptsDown };
    }

    private static final String FILE_EXTENSION = ".imu.neo_epitope.tsv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<NeoEpitopeFile> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<NeoEpitopeFile> neos) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(neos));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<NeoEpitopeFile> neos)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        neos.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<NeoEpitopeFile> fromLines(@NotNull List<String> lines)
    {
        if(lines.isEmpty())
            return Lists.newArrayList();

        final String header = lines.get(0);

        boolean skipSampleId = header.startsWith("SampleId");
        lines.remove(0);

        return lines.stream().map(x -> NeoEpitopeFile.fromString(x, skipSampleId)).collect(toList());
    }

    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("VariantType")
                .add("VariantInfo")
                .add("JunctionCopyNumber")
                .add("GeneIdUp")
                .add("GeneIdDown")
                .add("GeneNameUp")
                .add("GeneNameDown")
                .add("UpstreamAA")
                .add("DownstreamAA")
                .add("NovelAA")
                .add("NmdMin")
                .add("NmdMax")
                .add("UpTranscripts")
                .add("DownTranscripts")
                .add("CodingBaseUpPos")
                .add("CodingBaseDownPos")
                .add("CodingBasesUp")
                .add("CodingBasesDown")
                .add("WildtypeAA")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final NeoEpitopeFile neo)
    {
        StringJoiner sj = new StringJoiner(DELIMITER);

        sj.add(neo.VariantType.toString());
        sj.add(neo.VariantInfo);
        sj.add(String.format("%.4f", neo.CopyNumber));
        sj.add(neo.GeneIds[FS_UP]);
        sj.add(neo.GeneIds[FS_DOWN]);
        sj.add(neo.GeneNames[FS_UP]);
        sj.add(neo.GeneNames[FS_DOWN]);
        sj.add(neo.UpstreamAA);
        sj.add(neo.DownstreamAA);
        sj.add(neo.NovelAA);
        sj.add(String.valueOf(neo.NmdBases[0]));
        sj.add(String.valueOf(neo.NmdBases[1]));
        sj.add(neo.Transcripts[FS_UP]);
        sj.add(neo.Transcripts[FS_DOWN]);
        sj.add(String.valueOf(neo.CodingBasePositions[FS_UP]));
        sj.add(String.valueOf(neo.CodingBasePositions[FS_DOWN]));
        sj.add(neo.CodingBases[FS_UP]);
        sj.add(neo.CodingBases[FS_DOWN]);
        sj.add(neo.WildtypeAA);

        return sj.toString();
    }

    @NotNull
    public static NeoEpitopeFile fromString(@NotNull final String data, boolean skipSampleId)
    {
        final String[] values = data.split(DELIMITER, -1);

        int index = 0;

        if(skipSampleId)
            ++index;

        return new NeoEpitopeFile(
                NeoEpitopeType.valueOf(values[index++]), values[index++], Double.parseDouble(values[index++]),
                values[index++], values[index++], values[index++], values[index++],
                values[index++], values[index++], values[index++], values[index++],
                Integer.parseInt(values[index++]), Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]), Integer.parseInt(values[index++]),
                values[index++], values[index++], values[index++], values[index++]);
    }

    public static final String ITEM_DELIM = ";";
    public static final String VAR_INFO_DELIM = ":";
    public static final String FUSION_INFO_DELIM = "-";

    public static String pointMutationInfo(final String chromosome, int position, final String ref, final String alt)
    {
        return String.format("%s:%d:%s:%s", chromosome,  position, ref, alt);
    }

    public static String fusionInfo(final String[] chromosomes, final int[] positions, final byte[] orientations)
    {
        return String.format("%s:%d:%d-%s:%d:%d",
                chromosomes[FS_UP], positions[FS_UP], orientations[FS_UP],
                chromosomes[FS_DOWN], positions[FS_DOWN], orientations[FS_DOWN]);
    }

    public void extractLocationData(final String[] chromosomes, final int[] positions, final byte[] orientations)
    {
        if(VariantInfo.contains(FUSION_INFO_DELIM))
        {
            final String[] streamData = VariantInfo.split(FUSION_INFO_DELIM);
            chromosomes[FS_UP] = streamData[FS_UP].split(VAR_INFO_DELIM)[0];
            chromosomes[FS_DOWN] = streamData[FS_DOWN].split(VAR_INFO_DELIM)[0];
            positions[FS_UP] = Integer.parseInt(streamData[FS_UP].split(VAR_INFO_DELIM)[1]);
            positions[FS_DOWN] = Integer.parseInt(streamData[FS_DOWN].split(VAR_INFO_DELIM)[1]);
            orientations[FS_UP] = Byte.parseByte(streamData[FS_UP].split(VAR_INFO_DELIM)[2]);
            orientations[FS_DOWN] = Byte.parseByte(streamData[FS_DOWN].split(VAR_INFO_DELIM)[2]);
        }
        else
        {
            final String[] varData = VariantInfo.split(VAR_INFO_DELIM);
            chromosomes[FS_UP] = chromosomes[FS_DOWN] = varData[0];
            positions[FS_UP] = positions[FS_DOWN] = Integer.parseInt(varData[1]);
        }
    }

    public void extractTranscriptNames(final List<String> upNames, final List<String> downNames)
    {
        Arrays.stream(Transcripts[FS_UP].split(ITEM_DELIM)).forEach(x -> upNames.add(x));
        Arrays.stream(Transcripts[FS_DOWN].split(ITEM_DELIM)).forEach(x -> downNames.add(x));
    }

}
