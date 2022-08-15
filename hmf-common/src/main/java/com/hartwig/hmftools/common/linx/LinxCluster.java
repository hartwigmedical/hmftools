package com.hartwig.hmftools.common.linx;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.ImmutableLinxCluster;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxCluster
{
    public abstract int clusterId();
    public abstract String category();
    public abstract boolean synthetic();
    public abstract String resolvedType();
    public abstract int clusterCount();
    public abstract String clusterDesc();

    public static final String DELIMITER = "\t";

    private static final String FILE_EXTENSION = ".linx.clusters.tsv";
    private static final String GERMLINE_FILE_EXTENSION = ".linx.germline.clusters.tsv";

    public static String generateFilename(final String basePath, final String sample)
    {
        return generateFilename(basePath, sample, false);
    }

    public static String generateFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? GERMLINE_FILE_EXTENSION : FILE_EXTENSION);
    }

    public static List<LinxCluster> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<LinxCluster> clusters) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(clusters));
    }

    static List<String> toLines(final List<LinxCluster> clusters)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        clusters.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    static List<LinxCluster> fromLines(final List<String> lines)
    {
        final String header = lines.get(0);
        lines.remove(0);

        final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header,DELIMITER);

        List<LinxCluster> clusters = Lists.newArrayList();

        for(int i = 0; i < lines.size(); ++i)
        {
            String[] values = lines.get(i).split(DELIMITER);

            clusters.add(ImmutableLinxCluster.builder()
                    .clusterId(Integer.parseInt(values[fieldsIndexMap.get("clusterId")]))
                    .category(values[fieldsIndexMap.get("category")])
                    .synthetic(Boolean.parseBoolean(values[fieldsIndexMap.get("synthetic")]))
                    .resolvedType(values[fieldsIndexMap.get("resolvedType")])
                    .clusterCount(Integer.parseInt(values[fieldsIndexMap.get("clusterCount")]))
                    .clusterDesc(values[fieldsIndexMap.get("clusterDesc")])
                    .build());
        }

        return clusters;
    }

    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("clusterId")
                .add("category")
                .add("synthetic")
                .add("resolvedType")
                .add("clusterCount")
                .add("clusterDesc")
                .toString();
    }

    private static String toString(final LinxCluster cluster)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(cluster.clusterId()))
                .add(String.valueOf(cluster.category()))
                .add(String.valueOf(cluster.synthetic()))
                .add(String.valueOf(cluster.resolvedType()))
                .add(String.valueOf(cluster.clusterCount()))
                .add(String.valueOf(cluster.clusterDesc()))
                .toString();
    }
}
