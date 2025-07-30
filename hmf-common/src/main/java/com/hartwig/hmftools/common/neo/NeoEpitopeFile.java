package com.hartwig.hmftools.common.neo;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_PAIR;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public class NeoEpitopeFile
{
    public final int Id;
    public final NeoEpitopeType VariantType;
    public final String VariantInfo;
    public final double VariantCopyNumber;
    public final double CopyNumber;
    public final double SubclonalLikelihood;
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

    public static final String NEO_FILE_ID = "neo";

    public static final String FLD_NE_ID = "NeId";
    public static final String FLD_NE_VAR_TYPE = "VariantType";
    public static final String FLD_NE_VAR_INFO = "VariantInfo";
    public static final String FLD_NE_GENE_ID_UP = "GeneIdUp";
    public static final String FLD_NE_GENE_ID_DOWN = "GeneIdDown";
    public static final String FLD_NE_GENE_NAME_UP = "GeneNameUp";
    public static final String FLD_NE_GENE_NAME_DOWN = "GeneNameDown";
    public static final String FLD_NE_AA_UP = "UpstreamAA";
    public static final String FLD_NE_AA_DOWN = "DownstreamAA";
    public static final String FLD_NE_AA_NOVEL = "NovelAA";
    public static final String FLD_NE_TRANS_UP = "UpTranscripts";
    public static final String FLD_NE_TRANS_DOWN = "DownTranscripts";
    public static final String FLD_NE_NMD_MIN = "NmdMin";
    public static final String FLD_NE_NMD_MAX = "NmdMax";
    public static final String FLD_NE_CB_LEN_MIN = "CodingBasesLengthMin";
    public static final String FLD_NE_CB_LEN_MAX = "CodingBasesLengthMax";
    public static final String FLD_NE_FUSED_LEN = "FusedIntronLength";
    public static final String FLD_NE_SKIP_DONORS = "SkippedDonors";
    public static final String FLD_NE_SKIP_ACCEPTORS = "SkippedAcceptors";
    public static final String FLD_NE_VAR_CN = "VariantCopyNumber";
    public static final String FLD_NE_CN = "CopyNumber";
    public static final String FLD_NE_SC_LIKELIHOOD = "SubclonalLikelihood";

    public static final String VAR_INFO_DELIM = ":";
    private static final String FUSION_INFO_DELIM = ";";

    public NeoEpitopeFile(
            int id, final NeoEpitopeType varType, final String varInfo,
            final double variantCopyNumber, final double copyNumber, final double subclonalLikelihood,
            final String geneIdUp, final String geneIdDown, final String geneNameUp, final String geneNameDown,
            final String chrUp, final String chrDown, byte orientUp, byte orientDown,
            final String upAA, final String downAAs, final String novelAAs,
            int nmdBasesMin, int nmdBasesMax, int codingBasesLengthMin, int codingBasesLengthMax, int fusedIntronLength, int skippedDonors, int skippedAcceptors,
            final String transcriptsUp, final String transcriptsDown, final String wildtypeAAs,
            int codingBaseUpPosStart, int codingBaseUpPosEnd, final String codingBasesUp, final String codingBaseCigarUp,
            int codingBaseDownPosStart, int codingBaseDownPosEnd, final String codingBasesDown, final String codingBaseCigarDown)
    {
        Id = id;
        VariantType = varType;
        VariantInfo = varInfo;
        VariantCopyNumber = variantCopyNumber;
        CopyNumber = copyNumber;
        SubclonalLikelihood = subclonalLikelihood;
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
    }

    private static final String FILE_EXTENSION = ".neo.neo_data.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + FILE_EXTENSION;
    }

    public static List<NeoEpitopeFile> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, final List<NeoEpitopeFile> neos) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(neos));
    }

    private static List<String> toLines(final List<NeoEpitopeFile> neos)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        neos.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<NeoEpitopeFile> fromLines(final List<String> lines)
    {
        List<NeoEpitopeFile> neoepitopes = Lists.newArrayList();

        if(lines.isEmpty())
            return neoepitopes;

        final String header = lines.get(0);

        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        int neIdIndex = fieldsIndexMap.get(FLD_NE_ID);
        int varTypeIndex = fieldsIndexMap.get(FLD_NE_VAR_TYPE);
        int varInfoIndex = fieldsIndexMap.get(FLD_NE_VAR_INFO);
        int geneIdUpIndex = fieldsIndexMap.get(FLD_NE_GENE_ID_UP);
        int geneIdDownIndex = fieldsIndexMap.get(FLD_NE_GENE_ID_DOWN);
        int geneNameUpIndex = fieldsIndexMap.get(FLD_NE_GENE_NAME_UP);
        int geneNameDownIndex = fieldsIndexMap.get(FLD_NE_GENE_NAME_DOWN);
        int chrUpIndex = fieldsIndexMap.get("ChrUp");
        int chrDownIndex = fieldsIndexMap.get("ChrDown");
        int orientUpIndex = fieldsIndexMap.get("OrientUp");
        int orientDownIndex = fieldsIndexMap.get("OrientDown");
        int upAaIndex = fieldsIndexMap.get(FLD_NE_AA_UP);
        int downAaIndex = fieldsIndexMap.get(FLD_NE_AA_DOWN);
        int novelAaIndex = fieldsIndexMap.get(FLD_NE_AA_NOVEL);
        int nmdMinIndex = fieldsIndexMap.get(FLD_NE_NMD_MIN);
        int nmdMaxIndex = fieldsIndexMap.get(FLD_NE_NMD_MAX);
        int vcnIndex = fieldsIndexMap.get(FLD_NE_VAR_CN);
        int cnIndex = fieldsIndexMap.get(FLD_NE_CN);
        int sclIndex = fieldsIndexMap.get(FLD_NE_SC_LIKELIHOOD);
        int cbLenMinIndex = fieldsIndexMap.get(FLD_NE_CB_LEN_MIN);
        int cbLenMaxIndex = fieldsIndexMap.get(FLD_NE_CB_LEN_MAX);
        int feLenIndex = fieldsIndexMap.get(FLD_NE_FUSED_LEN);
        int skipDonIndex = fieldsIndexMap.get(FLD_NE_SKIP_DONORS);
        int skipAccIndex = fieldsIndexMap.get(FLD_NE_SKIP_ACCEPTORS);
        int transUpIndex = fieldsIndexMap.get(FLD_NE_TRANS_UP);
        int transDownIndex = fieldsIndexMap.get(FLD_NE_TRANS_DOWN);
        int wtAaIndex = fieldsIndexMap.get("WildtypeAA");
        int cbUpPosStartIndex = fieldsIndexMap.get("CodingBaseUpPosStart");
        int cbUpPosEndIndex = fieldsIndexMap.get("CodingBaseUpPosEnd");
        int cbUpIndex = fieldsIndexMap.get("CodingBasesUp");
        int cbCigUpIndex = fieldsIndexMap.get("CodingBaseCigarUp");
        int cbDownPosStartIndex = fieldsIndexMap.get("CodingBaseDownPosStart");
        int cbDownPosEndIndex = fieldsIndexMap.get("CodingBaseDownPosEnd");
        int cbDownIndex = fieldsIndexMap.get("CodingBasesDown");
        int cbCigDownIndex = fieldsIndexMap.get("CodingBaseCigarDown");

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            neoepitopes.add(new NeoEpitopeFile(
                    Integer.parseInt(values[neIdIndex]), NeoEpitopeType.valueOf(values[varTypeIndex]), values[varInfoIndex],
                    Double.parseDouble(values[vcnIndex]), Double.parseDouble(values[cnIndex]), Double.parseDouble(values[sclIndex]),
                    values[geneIdUpIndex], values[geneIdDownIndex], values[geneNameUpIndex], values[geneNameDownIndex],
                    values[chrUpIndex], values[chrDownIndex], Byte.parseByte(values[orientUpIndex]), Byte.parseByte(values[orientDownIndex]),
                    values[upAaIndex], values[downAaIndex], values[novelAaIndex], Integer.parseInt(values[nmdMinIndex]), Integer.parseInt(values[nmdMaxIndex]),
                    Integer.parseInt(values[cbLenMinIndex]), Integer.parseInt(values[cbLenMaxIndex]),
                    Integer.parseInt(values[feLenIndex]), Integer.parseInt(values[skipDonIndex]), Integer.parseInt(values[skipAccIndex]),
                    values[transUpIndex], values[transDownIndex], values[wtAaIndex],
                    Integer.parseInt(values[cbUpPosStartIndex]), Integer.parseInt(values[cbUpPosEndIndex]), values[cbUpIndex], values[cbCigUpIndex],
                    Integer.parseInt(values[cbDownPosStartIndex]), Integer.parseInt(values[cbDownPosEndIndex]), values[cbDownIndex], values[cbCigDownIndex]));
        }

        return neoepitopes;
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add(FLD_NE_ID)
                .add(FLD_NE_VAR_TYPE)
                .add(FLD_NE_VAR_INFO)
                .add(FLD_NE_VAR_CN)
                .add(FLD_NE_CN)
                .add(FLD_NE_SC_LIKELIHOOD)
                .add(FLD_NE_GENE_ID_UP)
                .add(FLD_NE_GENE_ID_DOWN)
                .add(FLD_NE_GENE_NAME_UP)
                .add(FLD_NE_GENE_NAME_DOWN)
                .add("ChrUp")
                .add("ChrDown")
                .add("OrientUp")
                .add("OrientDown")
                .add(FLD_NE_AA_UP)
                .add(FLD_NE_AA_DOWN)
                .add(FLD_NE_AA_NOVEL)
                .add(FLD_NE_NMD_MIN)
                .add(FLD_NE_NMD_MAX)
                .add(FLD_NE_CB_LEN_MIN)
                .add(FLD_NE_CB_LEN_MAX)
                .add(FLD_NE_FUSED_LEN)
                .add(FLD_NE_SKIP_DONORS)
                .add(FLD_NE_SKIP_ACCEPTORS)
                .add(FLD_NE_TRANS_UP)
                .add(FLD_NE_TRANS_DOWN)
                .add("WildtypeAA")
                .add("CodingBaseUpPosStart")
                .add("CodingBaseUpPosEnd")
                .add("CodingBasesUp")
                .add("CodingBaseCigarUp")
                .add("CodingBaseDownPosStart")
                .add("CodingBaseDownPosEnd")
                .add("CodingBasesDown")
                .add("CodingBaseCigarDown")
                .toString();
    }

    public static String toString(final NeoEpitopeFile neo)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(String.valueOf(neo.Id));
        sj.add(neo.VariantType.toString());
        sj.add(neo.VariantInfo);
        sj.add(String.format("%.4f", neo.VariantCopyNumber));
        sj.add(String.format("%.4f", neo.CopyNumber));
        sj.add(String.format("%.4f", neo.SubclonalLikelihood));
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

        return sj.toString();
    }

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
