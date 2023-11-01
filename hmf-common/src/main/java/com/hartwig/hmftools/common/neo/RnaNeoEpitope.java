package com.hartwig.hmftools.common.neo;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_ID;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_VAR_INFO;
import static com.hartwig.hmftools.common.neo.NeoEpitopeFile.FLD_NE_VAR_TYPE;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
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

import org.jetbrains.annotations.NotNull;

public class RnaNeoEpitope
{
    public final int Id;
    public final NeoEpitopeType VariantType;
    public final String VariantInfo;
    public final int FragmentCount;
    public final int[] BaseDepth;

    public static final String FLD_FRAG_COUNT = "FragmentCount";
    public static final String FLD_BASE_DEPTH_UP = "BaseDepthUp";
    public static final String FLD_BASE_DEPTH_DOWN = "BaseDepthDown";

    public RnaNeoEpitope(
            int id, final NeoEpitopeType varType, final String varInfo, int fragmentCount, int baseDepthUp, int baseDepthDown)
    {
        Id = id;
        VariantType = varType;
        VariantInfo = varInfo;
        FragmentCount = fragmentCount;
        BaseDepth = new int[] { baseDepthUp, baseDepthDown };
    }

    private static final String FILE_EXTENSION = ISF_FILE_ID + "neoepitope.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + FILE_EXTENSION;
    }

    public static List<RnaNeoEpitope> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<RnaNeoEpitope> neos) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(neos));
    }

    private static List<String> toLines(final List<RnaNeoEpitope> neos)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        neos.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    @NotNull
    private static List<RnaNeoEpitope> fromLines(List<String> lines)
    {
        if(lines.isEmpty())
            return Lists.newArrayList();

        final String header = lines.get(0);
        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        int neIdIndex = fieldsIndexMap.get(FLD_NE_ID);
        int varTypeIndex = fieldsIndexMap.get(FLD_NE_VAR_TYPE);
        int varInfoIndex = fieldsIndexMap.get(FLD_NE_VAR_INFO);
        int fragIndex = fieldsIndexMap.get(FLD_FRAG_COUNT);
        int baseDepthUpIndex = fieldsIndexMap.get(FLD_BASE_DEPTH_UP);
        int baseDepthDownIndex = fieldsIndexMap.get(FLD_BASE_DEPTH_DOWN);

        return lines.stream()
                .map(x -> RnaNeoEpitope.fromString(x, neIdIndex, varTypeIndex, varInfoIndex, fragIndex, baseDepthUpIndex, baseDepthDownIndex))
                .collect(toList());
    }

    public static String header()
    {
        return new StringJoiner(TSV_DELIM)
                .add(FLD_NE_ID)
                .add(FLD_NE_VAR_TYPE)
                .add(FLD_NE_VAR_INFO)
                .add(FLD_FRAG_COUNT)
                .add(FLD_BASE_DEPTH_UP)
                .add(FLD_BASE_DEPTH_DOWN)
                .toString();
    }

    @NotNull
    public static String toString(final RnaNeoEpitope neo)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(String.valueOf(neo.Id));
        sj.add(neo.VariantType.toString());
        sj.add(neo.VariantInfo);
        sj.add(String.valueOf(neo.FragmentCount));
        sj.add(String.valueOf(neo.BaseDepth[FS_UP]));
        sj.add(String.valueOf(neo.BaseDepth[FS_DOWN]));
        return sj.toString();
    }

    @NotNull
    public static RnaNeoEpitope fromString(
            final String data, int neIdIndex, int varTypeIndex, int varInfoIndex,
            int fragIndex, int baseDepthUpIndex, int baseDepthDownIndex)
    {
        final String[] values = data.split(TSV_DELIM, -1);

        return new RnaNeoEpitope(
                Integer.parseInt(values[neIdIndex]), NeoEpitopeType.valueOf(values[varTypeIndex]), values[varInfoIndex],
                Integer.parseInt(values[fragIndex]), Integer.parseInt(values[baseDepthUpIndex]), Integer.parseInt(values[baseDepthDownIndex]));
    }

    @NotNull
    public static RnaNeoEpitope fromString(final String data)
    {
        final String[] values = data.split(TSV_DELIM, -1);

        int index = 0;

        return new RnaNeoEpitope(
                Integer.parseInt(values[index++]), NeoEpitopeType.valueOf(values[index++]), values[index++],
                Integer.parseInt(values[index++]), Integer.parseInt(values[index++]), Integer.parseInt(values[index++]));
    }

}
