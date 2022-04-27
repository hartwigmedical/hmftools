package com.hartwig.hmftools.linx.visualiser.file;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber.DELIMITER;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.linx.types.ResolvedType;

import org.jetbrains.annotations.NotNull;

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

    @NotNull
    public static String generateFilename(@NotNull final String basePath, @NotNull final String sample)
    {
        return basePath + File.separator + sample + FILE_EXTENSION;
    }

    @NotNull
    public static List<VisSvData> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(@NotNull final String filename, @NotNull List<VisSvData> svDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svDataList));
    }

    @NotNull
    static List<String> toLines(@NotNull final List<VisSvData> cnDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        cnDataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    static List<VisSvData> fromLines(@NotNull List<String> lines)
    {
        return lines.stream().filter(x -> !x.startsWith("SampleId")).map(VisSvData::fromString).collect(toList());
    }

    @NotNull
    public static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("SampleId")
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

    @NotNull
    public static String toString(@NotNull final VisSvData svData)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(svData.SampleId))
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

    @NotNull
    private static VisSvData fromString(@NotNull final String tiData)
    {
        String[] values = tiData.split(DELIMITER);

        int index = 0;

        return new VisSvData(
                values[index++],
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                StructuralVariantType.valueOf(values[index++]),
                com.hartwig.hmftools.linx.types.ResolvedType.valueOf(values[index++]),
                Boolean.parseBoolean(values[index++]),
                values[index++],
                values[index++],
                Integer.parseInt(values[index++]),
                Integer.parseInt(values[index++]),
                Byte.parseByte(values[index++]),
                Byte.parseByte(values[index++]),
                values[index++],
                values[index++],
                Double.parseDouble(values[index++]),
                index < values.length ? Boolean.parseBoolean(values[index++]) : false);
    }

}
