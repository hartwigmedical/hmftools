package com.hartwig.hmftools.linx.visualiser.data;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class Fusions
{
    private static final String SEPARATOR = ",";

    @NotNull
    public static List<Fusion> fromFile(@NotNull final String fileName) throws IOException
    {
        return fromLines(Files.readAllLines(new File(fileName).toPath()));
    }

    @NotNull
    private static List<Fusion> fromLines(@NotNull final List<String> lines)
    {
        final Map<String, Integer> map = columnMap(lines.get(0));

        final List<Fusion> result = Lists.newArrayList();
        for (int i = 1; i < lines.size(); i++)
        {
            Fusion fusion = fromLine(map, lines.get(i));
            if (fusion.reportable())
            {
                result.add(fusion);
            }
        }

        return result;
    }

    @NotNull
    private static Map<String, Integer> columnMap(@NotNull final String header)
    {
        final Map<String, Integer> result = Maps.newHashMap();
        final String[] values = header.split(SEPARATOR);

        for (int i = 0; i < values.length; i++)
        {
            result.put(values[i], i);
        }

        return result;

    }

    private static Fusion fromLine(@NotNull final Map<String, Integer> map, @NotNull final String line)
    {

        String[] values = line.split(SEPARATOR);
        return ImmutableFusion.builder()
                .sampleId(values[map.get("SampleId")])
                .reportable(Boolean.valueOf(values[map.get("Reportable")]))
                .clusterId(Integer.valueOf(values[map.get("ClusterId")]))
                .chromosomeUp(values[map.get("ChrUp")])
                .positionUp(Long.valueOf(values[map.get("PosUp")]))
                .geneUp(values[map.get("GeneNameUp")])
                .strandUp(Integer.valueOf(values[map.get("StrandUp")]))
                .regionTypeUp(values[map.get("RegionTypeUp")])
                .exonUp(Integer.valueOf(values[map.get("ExonUp")]))
                .exonsSkippedUp(Integer.valueOf(values[map.get("ExonsSkippedUp")]))
                .chromosomeDown(values[map.get("ChrDown")])
                .positionDown(Long.valueOf(values[map.get("PosDown")]))
                .geneDown(values[map.get("GeneNameDown")])
                .strandDown(Integer.valueOf(values[map.get("StrandDown")]))
                .regionTypeDown(values[map.get("RegionTypeDown")])
                .exonDown(Integer.valueOf(values[map.get("ExonDown")]))
                .exonsSkippedDown(Integer.valueOf(values[map.get("ExonsSkippedDown")]))
                .build();
    }
}
