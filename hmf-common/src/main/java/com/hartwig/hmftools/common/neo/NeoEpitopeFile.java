package com.hartwig.hmftools.common.neo;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class NeoEpitopeFile
{
    public final int Id;
    public final NeoEpitopeType VariantType;
    public final String VariantInfo;
    public final double CopyNumber;
    public final String[] GeneIds;
    public final String[] GeneNames;

    public final String[] Chromosomes;
    public final byte[] Orientations;
    public final int[] NmdBases;
    public final int[] CodingBasesLength;
    public final String UpstreamAA;
    public final String DownstreamAA;
    public final String NovelAA;
    public final String[] Transcripts;
    public final String WildtypeAA;
    public final int FusedIntronLength;
    public final int[] SkippedAcceptorsDonors;
    public final int[][] CodingBasePositions;
    public final String[] CodingBases;
    public final String[] CodingBaseCigars;
    public final double[] CancerTpmTotal;
    public final double[] CohortTpmTotal;

    public static final String FLD_NE_ID = "NeId";
    public static final String FLD_NE_VAR_TYPE = "VariantType";
    public static final String FLD_NE_VAR_INFO = "VariantInfo";

    public static final String DELIMITER = ",";
    public static final String ITEM_DELIM = ";";
    public static final String VAR_INFO_DELIM = ":";
    private static final String FUSION_INFO_DELIM = ";";

    protected static final Logger NEO_LOGGER = LogManager.getLogger(NeoEpitopeFile.class);

    public NeoEpitopeFile(
            int id, final NeoEpitopeType varType, final String varInfo, final double copyNumber,
            final String geneIdUp, final String geneIdDown, final String geneNameUp, final String geneNameDown,
            final String chrUp, final String chrDown, byte orientUp, byte orientDown,
            final String upAA, final String downAAs, final String novelAAs,
            int nmdBasesMin, int nmdBasesMax, int codingBasesLengthMin, int codingBasesLengthMax, int fusedIntronLength, int skippedDonors, int skippedAcceptors,
            final String transcriptsUp, final String transcriptsDown, final String wildtypeAAs,
            int codingBaseUpPosStart, int codingBaseUpPosEnd, final String codingBasesUp, final String codingBaseCigarUp,
            int codingBaseDownPosStart, int codingBaseDownPosEnd, final String codingBasesDown, final String codingBaseCigarDown,
            double tmpCancerUp, double tmpCohortUp, double tmpCancerDown, double tmpCohortDown)
    {
        Id = id;
        VariantType = varType;
        VariantInfo = varInfo;
        CopyNumber = copyNumber;
        GeneIds = new String[] { geneIdUp, geneIdDown };
        GeneNames = new String[] { geneNameUp, geneNameDown };
        Chromosomes = new String[] { chrUp, chrDown };
        Orientations = new byte[] { orientUp, orientDown };
        NmdBases = new int[] { nmdBasesMin, nmdBasesMax };
        CodingBasesLength = new int[] { codingBasesLengthMin, codingBasesLengthMax };
        UpstreamAA = upAA;
        DownstreamAA = downAAs;
        NovelAA = novelAAs;
        WildtypeAA = wildtypeAAs;
        FusedIntronLength = fusedIntronLength;
        SkippedAcceptorsDonors = new int[] { skippedDonors, skippedAcceptors };
        Transcripts = new String[] { transcriptsUp, transcriptsDown };

        CodingBases = new String[] { codingBasesUp, codingBasesDown };
        CodingBaseCigars = new String[] { codingBaseCigarUp, codingBaseCigarDown };
        CodingBasePositions = new int[FS_PAIR][SE_PAIR];
        CodingBasePositions[FS_UP] = new int[] {codingBaseUpPosStart, codingBaseUpPosEnd};
        CodingBasePositions[FS_DOWN] = new int[] {codingBaseDownPosStart, codingBaseDownPosEnd};
        CancerTpmTotal = new double[] { tmpCancerUp, tmpCohortUp };
        CohortTpmTotal = new double[] { tmpCancerDown, tmpCohortDown };
    }

    private static final String FILE_EXTENSION = ".neo.neoepitopes.csv";
    private static final String OLD_FILE_EXTENSION = ".imu.neo_epitopes.csv";

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        String filename = basePath + File.separator + sample + OLD_FILE_EXTENSION;

        if(Files.exists(Paths.get(filename)))
            return filename;

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
        List<NeoEpitopeFile> neoepitopes = Lists.newArrayList();

        if(lines.isEmpty())
            return neoepitopes;

        final String header = lines.get(0);

        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);

        Integer sampleIndex = fieldsIndexMap.get(FLD_NE_ID);
        int neIdIndex = fieldsIndexMap.get(FLD_NE_ID);
        int varTypeIndex = fieldsIndexMap.get(FLD_NE_VAR_TYPE);
        int varInfoIndex = fieldsIndexMap.get(FLD_NE_VAR_INFO);
        int geneIdUpIndex = fieldsIndexMap.get("GeneIdUp");
        int geneIdDowwIndex = fieldsIndexMap.get("GeneIdDown");
        int geneNameUpIndex = fieldsIndexMap.get("GeneNameUp");
        int geneNameDownIndex = fieldsIndexMap.get("GeneNameDown");
        int chrUpIndex = fieldsIndexMap.get("ChrUp");
        int chrDownIndex = fieldsIndexMap.get("ChrDown");
        int orientUpIndex = fieldsIndexMap.get("OrientUp");
        int orientDownIndex = fieldsIndexMap.get("OrientDown");
        int upAaIndex = fieldsIndexMap.get("UpstreamAA");
        int downAaIndex = fieldsIndexMap.get("DownstreamAA");
        int novelAaIndex = fieldsIndexMap.get("NovelAA");
        int nmdMinIndex = fieldsIndexMap.get("NmdMin");
        int nmdMaxIndex = fieldsIndexMap.get("NmdMax");
        int jcnIndex = fieldsIndexMap.get("JunctionCopyNumber");
        int cbLenMinIndex = fieldsIndexMap.get("CodingBasesLengthMin");
        int cbLenMaxIndex = fieldsIndexMap.get("CodingBasesLengthMax");
        int feLenIndex = fieldsIndexMap.get("FusedIntronLength");
        int skipDonIndex = fieldsIndexMap.get("SkippedDonors");
        int skipAccIndex = fieldsIndexMap.get("SkippedAcceptors");
        int upTransIndex = fieldsIndexMap.get("UpTranscripts");
        int downTransIndex = fieldsIndexMap.get("DownTranscripts");
        int wtAaIndex = fieldsIndexMap.get("WildtypeAA");
        int cbUpPosStartIndex = fieldsIndexMap.get("CodingBaseUpPosStart");
        int cbUpPosEndIndex = fieldsIndexMap.get("CodingBaseUpPosEnd");
        int cbUpIndex = fieldsIndexMap.get("CodingBasesUp");
        int cbCigUpIndex = fieldsIndexMap.get("CodingBaseCigarUp");
        int cbDownPosStartIndex = fieldsIndexMap.get("CodingBaseDownPosStart");
        int cbDownPosEndIndex = fieldsIndexMap.get("CodingBaseDownPosEnd");
        int cbDownIndex = fieldsIndexMap.get("CodingBasesDown");
        int cbCigDownIndex = fieldsIndexMap.get("CodingBaseCigarDown");
        int tpmCanUpIndex = fieldsIndexMap.get("TpmCancerUp");
        int tpmCohUpIndex = fieldsIndexMap.get("TpmCohortUp");
        int tpmCanDownIndex = fieldsIndexMap.get("TpmCancerDown");
        int tpmCohDownIndex = fieldsIndexMap.get("TpmCohortDown");

        for(String line : lines)
        {
            String[] values = line.split(DELIMITER, -1);

            try
            {
                neoepitopes.add(new NeoEpitopeFile(
                        Integer.parseInt(values[neIdIndex]), NeoEpitopeType.valueOf(values[varTypeIndex]), values[varInfoIndex], Double.parseDouble(values[jcnIndex]),
                        values[geneIdUpIndex], values[geneIdDowwIndex], values[geneNameUpIndex], values[geneNameDownIndex],
                        values[chrUpIndex], values[chrDownIndex], Byte.parseByte(values[orientUpIndex]), Byte.parseByte(values[orientDownIndex]),
                        values[upAaIndex], values[downAaIndex], values[novelAaIndex], Integer.parseInt(values[nmdMinIndex]), Integer.parseInt(values[nmdMaxIndex]),
                        Integer.parseInt(values[cbLenMinIndex]), Integer.parseInt(values[cbLenMaxIndex]),
                        Integer.parseInt(values[feLenIndex]), Integer.parseInt(values[skipDonIndex]), Integer.parseInt(values[skipAccIndex]),
                        values[upTransIndex], values[downTransIndex], values[wtAaIndex],
                        Integer.parseInt(values[cbUpPosStartIndex]), Integer.parseInt(values[cbUpPosEndIndex]), values[cbUpIndex], values[cbCigUpIndex],
                        Integer.parseInt(values[cbDownPosStartIndex]), Integer.parseInt(values[cbDownPosEndIndex]), values[cbDownIndex], values[cbCigDownIndex],
                        Double.parseDouble(values[tpmCanUpIndex]), Double.parseDouble(values[tpmCohUpIndex]),
                        Double.parseDouble(values[tpmCanDownIndex]), Double.parseDouble(values[tpmCohDownIndex])));
            }
            catch(Exception e)
            {
                NEO_LOGGER.error("failed to read neoepitope line: {}", line);
            }

        }

        return neoepitopes;
    }

    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add(FLD_NE_ID)
                .add(FLD_NE_VAR_TYPE)
                .add(FLD_NE_VAR_INFO)
                .add("JunctionCopyNumber")
                .add("GeneIdUp")
                .add("GeneIdDown")
                .add("GeneNameUp")
                .add("GeneNameDown")
                .add("ChrUp")
                .add("ChrDown")
                .add("OrientUp")
                .add("OrientDown")
                .add("UpstreamAA")
                .add("DownstreamAA")
                .add("NovelAA")
                .add("NmdMin")
                .add("NmdMax")
                .add("CodingBasesLengthMin")
                .add("CodingBasesLengthMax")
                .add("FusedIntronLength")
                .add("SkippedDonors")
                .add("SkippedAcceptors")
                .add("UpTranscripts")
                .add("DownTranscripts")
                .add("WildtypeAA")
                .add("CodingBaseUpPosStart")
                .add("CodingBaseUpPosEnd")
                .add("CodingBasesUp")
                .add("CodingBaseCigarUp")
                .add("CodingBaseDownPosStart")
                .add("CodingBaseDownPosEnd")
                .add("CodingBasesDown")
                .add("CodingBaseCigarDown")
                .add("TpmCancerUp")
                .add("TpmCohortUp")
                .add("TpmCancerDown")
                .add("TpmCohortDown")
                .toString();
    }

    @NotNull
    public static String toString(@NotNull final NeoEpitopeFile neo)
    {
        StringJoiner sj = new StringJoiner(DELIMITER);

        sj.add(String.valueOf(neo.Id));
        sj.add(neo.VariantType.toString());
        sj.add(neo.VariantInfo);
        sj.add(String.format("%.4f", neo.CopyNumber));
        sj.add(neo.GeneIds[FS_UP]);
        sj.add(neo.GeneIds[FS_DOWN]);
        sj.add(neo.GeneNames[FS_UP]);
        sj.add(neo.GeneNames[FS_DOWN]);
        sj.add(neo.Chromosomes[FS_UP]);
        sj.add(neo.Chromosomes[FS_DOWN]);
        sj.add(String.valueOf(neo.Orientations[FS_UP]));
        sj.add(String.valueOf(neo.Orientations[FS_DOWN]));
        sj.add(neo.UpstreamAA);
        sj.add(neo.DownstreamAA);
        sj.add(neo.NovelAA);
        sj.add(String.valueOf(neo.NmdBases[0]));
        sj.add(String.valueOf(neo.NmdBases[1]));
        sj.add(String.valueOf(neo.CodingBasesLength[0]));
        sj.add(String.valueOf(neo.CodingBasesLength[1]));
        sj.add(String.valueOf(neo.FusedIntronLength));
        sj.add(String.valueOf(neo.SkippedAcceptorsDonors[FS_UP]));
        sj.add(String.valueOf(neo.SkippedAcceptorsDonors[FS_DOWN]));
        sj.add(neo.Transcripts[FS_UP]);
        sj.add(neo.Transcripts[FS_DOWN]);
        sj.add(neo.WildtypeAA);
        sj.add(String.valueOf(neo.CodingBasePositions[FS_UP][SE_START]));
        sj.add(String.valueOf(neo.CodingBasePositions[FS_UP][SE_END]));
        sj.add(neo.CodingBases[FS_UP]);
        sj.add(neo.CodingBaseCigars[FS_UP]);
        sj.add(String.valueOf(neo.CodingBasePositions[FS_DOWN][SE_START]));
        sj.add(String.valueOf(neo.CodingBasePositions[FS_DOWN][SE_END]));
        sj.add(neo.CodingBases[FS_DOWN]);
        sj.add(neo.CodingBaseCigars[FS_DOWN]);
        sj.add(String.format("%6.3e", neo.CancerTpmTotal[FS_UP]));
        sj.add(String.format("%6.3e", neo.CohortTpmTotal[FS_UP]));
        sj.add(String.format("%6.3e", neo.CancerTpmTotal[FS_DOWN]));
        sj.add(String.format("%6.3e", neo.CohortTpmTotal[FS_DOWN]));

        return sj.toString();
    }

    /*
    @NotNull
    public static NeoEpitopeFile fromString(@NotNull final String data, boolean skipSampleId)
    {
        final String[] values = data.split(DELIMITER, -1);

        int index = 0;

        if(skipSampleId)
            ++index;

        return new NeoEpitopeFile(
                Integer.parseInt(values[index++]), NeoEpitopeType.valueOf(values[index++]), values[index++], Double.parseDouble(values[index++]),
                values[index++], values[index++], values[index++], values[index++],
                values[index++], values[index++], Byte.parseByte(values[index++]), Byte.parseByte(values[index++]),
                values[index++], values[index++], values[index++], Integer.parseInt(values[index++]), Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]), Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]), Integer.parseInt(values[index++]), Integer.parseInt(values[index++]),
                values[index++], values[index++], values[index++],
                Integer.parseInt(values[index++]), Integer.parseInt(values[index++]), values[index++], values[index++],
                Integer.parseInt(values[index++]), Integer.parseInt(values[index++]), values[index++], values[index++],
                Double.parseDouble(values[index++]), Double.parseDouble(values[index++]),
                Double.parseDouble(values[index++]), Double.parseDouble(values[index++]));
    }
    */

    public static String pointMutationInfo(final String chromosome, int position, final String ref, final String alt)
    {
        return String.format("%s:%d:%s:%s", chromosome,  position, ref, alt);
    }

    public static String fusionInfo(final String[] chromosomes, final int[] positions, final byte[] orientations)
    {
        return String.format("%s:%d:%d;%s:%d:%d",
                chromosomes[FS_UP], positions[FS_UP], orientations[FS_UP],
                chromosomes[FS_DOWN], positions[FS_DOWN], orientations[FS_DOWN]);
    }

    public void extractLocationData(final String[] chromosomes, final int[] positions, final byte[] orientations)
    {
        extractLocationData(VariantInfo, chromosomes, positions, orientations);
    }

    public static void extractLocationData(final String variantInfo, final String[] chromosomes, final int[] positions, final byte[] orientations)
    {
        if(variantInfo.contains(FUSION_INFO_DELIM))
        {
            final String[] streamData = variantInfo.split(FUSION_INFO_DELIM);
            chromosomes[FS_UP] = streamData[FS_UP].split(VAR_INFO_DELIM)[0];
            chromosomes[FS_DOWN] = streamData[FS_DOWN].split(VAR_INFO_DELIM)[0];
            positions[FS_UP] = Integer.parseInt(streamData[FS_UP].split(VAR_INFO_DELIM)[1]);
            positions[FS_DOWN] = Integer.parseInt(streamData[FS_DOWN].split(VAR_INFO_DELIM)[1]);
            orientations[FS_UP] = Byte.parseByte(streamData[FS_UP].split(VAR_INFO_DELIM)[2]);
            orientations[FS_DOWN] = Byte.parseByte(streamData[FS_DOWN].split(VAR_INFO_DELIM)[2]);
        }
        else
        {
            final String[] varData = variantInfo.split(VAR_INFO_DELIM);
            chromosomes[FS_UP] = chromosomes[FS_DOWN] = varData[0];
            positions[FS_UP] = positions[FS_DOWN] = Integer.parseInt(varData[1]);
        }
    }

    public void extractTranscriptNames(final List<String> upNames, final List<String> downNames)
    {
        extractTranscriptNames(Transcripts[FS_UP], Transcripts[FS_DOWN], upNames, downNames);
    }

    public static void extractTranscriptNames(final String transUp, final String transDown, final List<String> upNames, final List<String> downNames)
    {
        Arrays.stream(transUp.split(ITEM_DELIM)).forEach(x -> upNames.add(x));
        Arrays.stream(transDown.split(ITEM_DELIM)).forEach(x -> downNames.add(x));
    }

}
