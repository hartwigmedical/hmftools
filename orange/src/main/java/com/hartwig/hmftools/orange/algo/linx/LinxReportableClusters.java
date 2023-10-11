package com.hartwig.hmftools.orange.algo.linx;

import static java.util.stream.Collectors.toList;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.getIntValue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxData;

import org.jetbrains.annotations.NotNull;

// TODO this is redundant with linx.visualizer.SampleData.findReportableClusters(),
//  remove when available via shared code.
public final class LinxReportableClusters
{
    @NotNull
    public static Set<Integer> findReportableClusters(@NotNull LinxData linxData, @NotNull String dataDir, @NotNull String sample)
            throws IOException
    {
        Set<Integer> clusterIds = Sets.newHashSet();

        clusterIds.addAll(loadFusionClusters(dataDir, sample));

        List<LinxBreakend> breakends = linxData.reportableSomaticBreakends();
        List<Integer> svIds = breakends.stream().map(LinxBreakend::svId).collect(toList());

        Map<Integer, Integer> svToCluster = loadSvToCluster(dataDir, sample);
        svIds.stream().filter(svToCluster::containsKey).map(svToCluster::get).forEach(clusterIds::add);

        linxData.somaticDrivers().stream().filter(x -> x.clusterId() >= 0).forEach(x -> clusterIds.add(x.clusterId()));
        return clusterIds;
    }

    private static final String VIS_FUSION_FILE_EXTENSION = ".linx.vis_fusion.tsv";

    @NotNull
    private static Set<Integer> loadFusionClusters(@NotNull String dataDir, @NotNull String sample) throws IOException
    {
        String filename = dataDir + File.separator + sample + VIS_FUSION_FILE_EXTENSION;
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        if(lines.isEmpty())
        {
            throw new IllegalStateException(String.format("File lacks header: %s", filename));
        }

        Set<Integer> clusterIds = Sets.newHashSet();
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);
            clusterIds.add(getIntValue(fieldsIndexMap, "ClusterId", 0, values));
        }

        return clusterIds;
    }

    private static final String VIS_SV_DATA_FILE_EXTENSION = ".linx.vis_sv_data.tsv";

    @NotNull
    private static Map<Integer, Integer> loadSvToCluster(@NotNull String dataDir, @NotNull String sample) throws IOException
    {
        String filename = dataDir + File.separator + sample + VIS_SV_DATA_FILE_EXTENSION;
        List<String> lines = Files.readAllLines(new File(filename).toPath());
        if(lines.isEmpty())
        {
            throw new IllegalStateException(String.format("File lacks header: %s", filename));
        }

        Map<Integer, Integer> svTocluster = Maps.newHashMap();
        String header = lines.get(0);
        Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        lines.remove(0);

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM);
            int clusterId = getIntValue(fieldsIndexMap, "ClusterId", 0, values);
            int svId = getIntValue(fieldsIndexMap, "SvId", 0, values);

            if(!svTocluster.containsKey(svId))
            {
                svTocluster.put(svId, clusterId);
            }
        }

        return svTocluster;
    }
}
