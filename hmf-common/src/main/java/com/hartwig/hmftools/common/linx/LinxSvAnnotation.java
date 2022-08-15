package com.hartwig.hmftools.common.linx;

import static com.hartwig.hmftools.common.linx.LinxCluster.DELIMITER;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.immutables.value.Value;

@Value.Immutable
public abstract class LinxSvAnnotation
{
    public abstract String vcfId();
    public abstract int svId();
    public abstract int clusterId();
    public abstract String clusterReason();
    public abstract boolean fragileSiteStart();
    public abstract boolean fragileSiteEnd();
    public abstract boolean isFoldback();
    public abstract String lineTypeStart();
    public abstract String lineTypeEnd();
    public abstract double junctionCopyNumberMin();
    public abstract double junctionCopyNumberMax();
    public abstract String geneStart();
    public abstract String geneEnd();
    public abstract int localTopologyIdStart();
    public abstract int localTopologyIdEnd();
    public abstract String localTopologyStart();
    public abstract String localTopologyEnd();
    public abstract int localTICountStart();
    public abstract int localTICountEnd();

    private static final String FILE_EXTENSION = ".linx.svs.tsv";
    private static final String GERMLINE_FILE_EXTENSION = ".linx.germline.svs.tsv";

    private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.0000");

    public static String generateFilename(final String basePath, final String sample)
    {
        return generateFilename(basePath, sample, false);
    }

    public static String generateFilename(final String basePath, final String sample, boolean isGermline)
    {
        return basePath + File.separator + sample + (isGermline ? GERMLINE_FILE_EXTENSION : FILE_EXTENSION);
    }

    public static List<LinxSvAnnotation> read(final String filePath) throws IOException
    {
        return fromLines(Files.readAllLines(new File(filePath).toPath()));
    }

    public static void write(final String filename, List<LinxSvAnnotation> svDataList) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(svDataList));
    }

    private static List<String> toLines(final List<LinxSvAnnotation> svDataList)
    {
        final List<String> lines = Lists.newArrayList();
        lines.add(header());
        svDataList.stream().map(x -> toString(x)).forEach(lines::add);
        return lines;
    }

    private static List<LinxSvAnnotation> fromLines(final List<String> lines)
    {
        String header = lines.get(0);
        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, DELIMITER);

        List<LinxSvAnnotation> annotations = Lists.newArrayList();

        for(int i = 1; i < lines.size(); ++i)
        {
            String[] values = lines.get(i).split(DELIMITER);

            annotations.add(ImmutableLinxSvAnnotation.builder()
                    .vcfId(values[fieldsIndexMap.get("vcfId")])
                    .svId(Integer.parseInt(values[fieldsIndexMap.get("svId")]))
                    .clusterId(Integer.parseInt(values[fieldsIndexMap.get("clusterId")]))
                    .clusterReason(values[fieldsIndexMap.get("clusterReason")])
                    .fragileSiteStart(Boolean.parseBoolean(values[fieldsIndexMap.get("fragileSiteStart")]))
                    .fragileSiteEnd(Boolean.parseBoolean(values[fieldsIndexMap.get("fragileSiteEnd")]))
                    .isFoldback(Boolean.parseBoolean(values[fieldsIndexMap.get("isFoldback")]))
                    .lineTypeStart(values[fieldsIndexMap.get("lineTypeStart")])
                    .lineTypeEnd(values[fieldsIndexMap.get("lineTypeEnd")])
                    .junctionCopyNumberMin(Double.parseDouble(values[fieldsIndexMap.get("junctionCopyNumberMin")]))
                    .junctionCopyNumberMax(Double.parseDouble(values[fieldsIndexMap.get("junctionCopyNumberMax")]))
                    .geneStart(values[fieldsIndexMap.get("geneStart")])
                    .geneEnd(values[fieldsIndexMap.get("geneEnd")])
                    .localTopologyIdStart(Integer.parseInt(values[fieldsIndexMap.get("localTopologyIdStart")]))
                    .localTopologyIdEnd(Integer.parseInt(values[fieldsIndexMap.get("localTopologyIdEnd")]))
                    .localTopologyStart(values[fieldsIndexMap.get("localTopologyStart")])
                    .localTopologyEnd(values[fieldsIndexMap.get("localTopologyEnd")])
                    .localTICountStart(Integer.parseInt(values[fieldsIndexMap.get("localTICountStart")]))
                    .localTICountEnd(Integer.parseInt(values[fieldsIndexMap.get("localTICountEnd")]))
                    .build());
        }

        return annotations;
    }

    private static String header()
    {
        return new StringJoiner(DELIMITER)
                .add("vcfId")
                .add("svId")
                .add("clusterId")
                .add("clusterReason")
                .add("fragileSiteStart")
                .add("fragileSiteEnd")
                .add("isFoldback")
                .add("lineTypeStart")
                .add("lineTypeEnd")
                .add("junctionCopyNumberMin")
                .add("junctionCopyNumberMax")
                .add("geneStart")
                .add("geneEnd")
                .add("localTopologyIdStart")
                .add("localTopologyIdEnd")
                .add("localTopologyStart")
                .add("localTopologyEnd")
                .add("localTICountStart")
                .add("localTICountEnd")
                .toString();
    }

    private static String toString(final LinxSvAnnotation svData)
    {
        return new StringJoiner(DELIMITER)
                .add(String.valueOf(svData.vcfId()))
                .add(String.valueOf(svData.svId()))
                .add(String.valueOf(svData.clusterId()))
                .add(String.valueOf(svData.clusterReason()))
                .add(String.valueOf(svData.fragileSiteStart()))
                .add(String.valueOf(svData.fragileSiteEnd()))
                .add(String.valueOf(svData.isFoldback()))
                .add(String.valueOf(svData.lineTypeStart()))
                .add(String.valueOf(svData.lineTypeEnd()))
                .add(DECIMAL_FORMAT.format(svData.junctionCopyNumberMin()))
                .add(DECIMAL_FORMAT.format(svData.junctionCopyNumberMax()))
                .add(String.valueOf(svData.geneStart()))
                .add(String.valueOf(svData.geneEnd()))
                .add(String.valueOf(svData.localTopologyIdStart()))
                .add(String.valueOf(svData.localTopologyIdEnd()))
                .add(String.valueOf(svData.localTopologyStart()))
                .add(String.valueOf(svData.localTopologyEnd()))
                .add(String.valueOf(svData.localTICountStart()))
                .add(String.valueOf(svData.localTICountEnd()))
                .toString();
    }
}
