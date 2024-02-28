package com.hartwig.hmftools.linx.visualiser.file;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getBoolValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getDoubleValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getValue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.types.ResolvedType;

public class VisSvData
{
    public final String SampleId;
    public final int ClusterId;
    public final int ChainId;
    public final int SvId;
    public final StructuralVariantType Type;
    public final ResolvedType ClusterResolvedType;
    public final boolean IsSynthetic;
    public final String ChrStart;
    public final String InfoStart;
    public final String ChrEnd;
    public final String InfoEnd;
    public int PosStart;
    public int PosEnd;
    public final byte OrientStart;
    public final byte OrientEnd;
    public final double JCN;
    public final boolean InDoubleMinute;

    public int Frame;

    public static final String INFO_TYPE_NORMAL = "NORMAL";
    public static final String INFO_TYPE_FOLDBACK = "FOLDBACK";

    public VisSvData(
            final String sampleId, int clusterId, int chainId, int svId, final StructuralVariantType type, final ResolvedType resolvedType,
            boolean isSynthetic, final String chrStart, final String chrEnd, int posStart, int posEnd, byte orientStart, byte orientEnd,
            final String infoStart, final String infoEnd, double jcn, boolean inDM)
    {
        SampleId = sampleId;
        ClusterId = clusterId;
        ChainId = chainId;
        ClusterResolvedType = resolvedType;
        IsSynthetic = isSynthetic;
        SvId = svId;
        Type = type;
        ChrStart = chrStart;
        ChrEnd = chrEnd;
        PosStart = posStart;
        PosEnd = posEnd;
        OrientStart = orientStart;
        OrientEnd = orientEnd;
        InfoStart = infoStart;
        InfoEnd = infoEnd;
        JCN = jcn;
        InDoubleMinute = inDM;
        Frame = 0;
    }

    public static VisSvData from(final VisSvData other)
    {
        VisSvData newData = new VisSvData(
                other.SampleId, other.ClusterId, other.ChainId, other.SvId, other.Type, other.ClusterResolvedType,
                other.IsSynthetic, other.ChrStart, other.ChrEnd, other.PosStart, other.PosEnd, other.OrientStart, other.OrientEnd,
                other.InfoStart, other.InfoEnd, other.JCN, other.InDoubleMinute);

        newData.Frame = other.Frame;

        return newData;
    }

    public boolean connectorsOnly(boolean showSimpleSvSegments)
    {
        return (!showSimpleSvSegments && isSimpleSV()) || isLineElement();
    }
    public boolean isSimpleSV()
    {
        return ClusterResolvedType.isSimple() && !IsSynthetic;
    }
    public boolean isLineElement()
    {
        return ClusterResolvedType == ResolvedType.LINE;
    }
    public boolean isValidStart()
    {
        return HumanChromosome.contains(ChrStart);
    }
    public boolean isValidEnd()
    {
        return HumanChromosome.contains(ChrEnd);
    }

    private static final String FILE_EXTENSION = ".linx.vis_sv_data.tsv";
    private static final String GERMLINE_FILE_EXTENSION = ".linx.germline.vis_sv_data.tsv";

    public static String generateFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? GERMLINE_FILE_EXTENSION : FILE_EXTENSION);
    }

    public static List<VisSvData> read(final String filePath) throws IOException
    {
        return fromLines(Files. readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<VisSvData> svDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svDataList));
    }

    private static List<String> toLines(final List<VisSvData> cnDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        cnDataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<VisSvData> fromLines(final List<String> lines)
    {
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        List<VisSvData> data = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);

            data.add(new VisSvData(
                    getValue(fieldsIndexMap, FLD_SAMPLE_ID, "", values),
                    getIntValue(fieldsIndexMap, "ClusterId", 0, values),
                    getIntValue(fieldsIndexMap, "ChainId", 0, values),
                    getIntValue(fieldsIndexMap, "SvId", 0, values),
                    StructuralVariantType.valueOf(values[fieldsIndexMap.get("Type")]),
                    ResolvedType.valueOf(values[fieldsIndexMap.get("ResolvedType")]),
                    getBoolValue(fieldsIndexMap, "IsSynthetic", false, values),
                    getValue(fieldsIndexMap, "ChrStart", "", values),
                    getValue(fieldsIndexMap, "ChrEnd", "", values),
                    getIntValue(fieldsIndexMap, "PosStart", 0, values),
                    getIntValue(fieldsIndexMap, "PosEnd", 0, values),
                    Byte.parseByte(values[fieldsIndexMap.get("OrientStart")]),
                    Byte.parseByte(values[fieldsIndexMap.get("OrientEnd")]),
                    getValue(fieldsIndexMap, "InfoStart", "", values),
                    getValue(fieldsIndexMap, "InfoEnd", "", values),
                    getDoubleValue(fieldsIndexMap, "JunctionCopyNumber", 0, values),
                    getBoolValue(fieldsIndexMap, "InDoubleMinute", false, values)));
        }

        return data;
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add("ClusterId")
                .add("ChainId")
                .add("SvId")
                .add("Type")
                .add("ResolvedType")
                .add("IsSynthetic")
                .add("ChrStart")
                .add("ChrEnd")
                .add("PosStart")
                .add("PosEnd")
                .add("OrientStart")
                .add("OrientEnd")
                .add("InfoStart")
                .add("InfoEnd")
                .add("JunctionCopyNumber")
                .add("InDoubleMinute")
                .toString();
    }

    public static String toString(final VisSvData svData)
    {
        return new StringJoiner(TSV_DELIM)
                .add(String.valueOf(svData.ClusterId))
                .add(String.valueOf(svData.ChainId))
                .add(String.valueOf(svData.SvId))
                .add(String.valueOf(svData.Type))
                .add(String.valueOf(svData.ClusterResolvedType))
                .add(String.valueOf(svData.IsSynthetic))
                .add(String.valueOf(svData.ChrStart))
                .add(String.valueOf(svData.ChrEnd))
                .add(String.valueOf(svData.PosStart))
                .add(String.valueOf(svData.PosEnd))
                .add(String.valueOf(svData.OrientStart))
                .add(String.valueOf(svData.OrientEnd))
                .add(String.valueOf(svData.InfoStart))
                .add(String.valueOf(svData.InfoEnd))
                .add(String.format("%.4f",svData.JCN))
                .add(String.valueOf(svData.InDoubleMinute))
                .toString();
    }
}
