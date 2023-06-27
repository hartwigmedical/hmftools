package com.hartwig.hmftools.common.neo;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public class NeoEpitopeFusion
{
    public final String[] GeneIds;
    public final String[] GeneNames;
    public final String[] Chromosomes;
    public final int[] Positions;
    public final byte[] Orientations;
    public final int[] SvIds;
    public final double JunctionCopyNumber;
    public final double CopyNumber;
    public final String InsertSequence;
    public final int ChainLength;

    public final String[] Transcripts;

    public NeoEpitopeFusion(
            final String geneIdUp, final String geneNameUp, final String chrUp, int posUp, byte orientUp, int svIdUp,
            final String geneIdDown, final String geneNameDown, final String chrDown, int posDown, byte orientDown, int svIdDown,
            final double junctionCopyNumber, final double copyNumber, final String insertSequence, final int chainLength, final String[] transcripts)
    {
        GeneIds = new String[] { geneIdUp, geneIdDown };
        GeneNames = new String[] { geneNameUp, geneNameDown };
        Chromosomes = new String[] { chrUp, chrDown };
        Positions = new int[] { posUp, posDown };
        Orientations = new byte[] { orientUp, orientDown };
        SvIds = new int[] { svIdUp, svIdDown };
        JunctionCopyNumber = junctionCopyNumber;
        CopyNumber = copyNumber;
        InsertSequence = insertSequence;
        ChainLength = chainLength;
        Transcripts = transcripts;
    }

    private static final String FILE_EXTENSION = ".linx.neoepitope.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + FILE_EXTENSION;
    }

    public static List<NeoEpitopeFusion> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<NeoEpitopeFusion> fusions) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(fusions));
    }

    private static List<String> toLines(final List<NeoEpitopeFusion> fusions)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        fusions.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<NeoEpitopeFusion> fromLines(final List<String> lines)
    {
        if(lines.isEmpty())
            return Lists.newArrayList();

        final String header = lines.get(0);
        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        List<NeoEpitopeFusion> fusions = Lists.newArrayList();

        int geneIdUpIndex = fieldsIndexMap.get("geneIdUp");
        int geneIdDownIndex = fieldsIndexMap.get("geneIdDown");
        int geneNameUpIndex = fieldsIndexMap.get("geneNameUp");
        int geneNameDownIndex = fieldsIndexMap.get("geneNameDown");
        int chrUpIndex = fieldsIndexMap.get("chromosomeUp");
        int chrDownIndex = fieldsIndexMap.get("chromosomeDown");
        int posUpIndex = fieldsIndexMap.get("positionUp");
        int posDownIndex = fieldsIndexMap.get("positionDown");
        int orientUpIndex = fieldsIndexMap.get("orientationUp");
        int orientDownIndex = fieldsIndexMap.get("orientationDown");
        int svIdUpIndex = fieldsIndexMap.get("svIdUp");
        int svIdDownIndex = fieldsIndexMap.get("svIdDown");
        int jcnIndex = fieldsIndexMap.get("junctionCopyNumber");
        int cnIndex = fieldsIndexMap.get("copyNumber");
        int insSeqIndex = fieldsIndexMap.get("insertSeq");
        int chainLengthIndex = fieldsIndexMap.get("chainLength");
        int transUpIndex = fieldsIndexMap.get("transcriptsUp");
        int transDownIndex = fieldsIndexMap.get("transcriptsDown");

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            fusions.add(new NeoEpitopeFusion(
                    values[geneIdUpIndex], values[geneNameUpIndex], values[chrUpIndex],
                    Integer.parseInt(values[posUpIndex]), Byte.parseByte(values[orientUpIndex]), Integer.parseInt(values[svIdUpIndex]),
                    values[geneIdDownIndex], values[geneNameDownIndex], values[chrDownIndex],
                    Integer.parseInt(values[posDownIndex]), Byte.parseByte(values[orientDownIndex]), Integer.parseInt(values[svIdDownIndex]),
                    Double.parseDouble(values[jcnIndex]), Double.parseDouble(values[cnIndex]), values[insSeqIndex],
                    Integer.parseInt(values[chainLengthIndex]), new String[] {values[transUpIndex], values[transDownIndex]}));
        }

        return fusions;
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
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
                .add("junctionCopyNumber")
                .add("copyNumber")
                .add("insertSeq")
                .add("chainLength")
                .add("transcriptsUp")
                .add("transcriptsDown")
                .toString();
    }

    public static String toString(final NeoEpitopeFusion fusion)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            sj.add(String.valueOf(fusion.GeneIds[fs]));
            sj.add(String.valueOf(fusion.GeneNames[fs]));
            sj.add(String.valueOf(fusion.Chromosomes[fs]));
            sj.add(String.valueOf(fusion.Positions[fs]));
            sj.add(String.valueOf(fusion.Orientations[fs]));
            sj.add(String.valueOf(fusion.SvIds[fs]));
        }

        sj.add(String.format("%.4f",fusion.JunctionCopyNumber));
        sj.add(String.format("%.4f",fusion.CopyNumber));
        sj.add(fusion.InsertSequence);
        sj.add(String.valueOf(fusion.ChainLength));
        sj.add(String.valueOf(fusion.Transcripts[FS_UP]));
        sj.add(String.valueOf(fusion.Transcripts[FS_DOWN]));

        return sj.toString();
    }
}
