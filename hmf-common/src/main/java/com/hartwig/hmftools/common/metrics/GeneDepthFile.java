package com.hartwig.hmftools.common.metrics;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.metrics.BamMetricsCommon.BAM_METRICS_FILE_ID;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public final class GeneDepthFile
{
    public static final String FLD_MV_LIKELIHOOD = "MissedVariantLikelihood";

    public static String generateGeneCoverageFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + BAM_METRICS_FILE_ID + ".gene_coverage.tsv";
    }

    public static String generateExonCoverageFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + BAM_METRICS_FILE_ID + ".exon_coverage.tsv";
    }

    public static void write(final String filename, final List<GeneDepth> depths, final List<Integer> depthBuckets) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(depths, depthBuckets));
    }

    public static class GeneDepthComparator implements Comparator<GeneDepth>
    {
        public int compare(final GeneDepth first, final GeneDepth second)
        {
            return first.Gene.compareTo(second.Gene);
        }
    }

    public static List<GeneDepth> read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String header = lines.get(0);
        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
        int geneIndex = fieldsIndexMap.get(FLD_GENE_NAME);
        int posStartIndex = fieldsIndexMap.get(FLD_POS_START);
        int posEndIndex = fieldsIndexMap.get(FLD_POS_END);
        int mvlIndex = fieldsIndexMap.get(FLD_MV_LIKELIHOOD);
        int depthStartIndex = mvlIndex + 1;

        List<GeneDepth> geneDepths = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            int depthColumns = values.length - depthStartIndex;

            int[] depthCounts = new int[depthColumns];

            for(int i = 0; i < depthColumns; ++i)
            {
                int colIndex = depthStartIndex + i;
                depthCounts[i] = Integer.parseInt(values[colIndex]);
            }

            GeneDepth geneDepth = new GeneDepth(
                    values[geneIndex], values[chrIndex], Integer.parseInt(values[posStartIndex]), Integer.parseInt(values[posEndIndex]),
                    Double.parseDouble(values[mvlIndex]), depthCounts);

            geneDepths.add(geneDepth);
        }

        return geneDepths;
    }

    private static List<String> toLines(final List<GeneDepth> depths, final List<Integer> depthBuckets)
    {
        if(!depths.isEmpty())
        {
            Collections.sort(depths, new GeneDepthComparator());
            final List<String> lines = Lists.newArrayList();
            lines.add(geneCoverageHeader(depthBuckets));
            depths.stream().map(GeneDepthFile::toString).forEach(lines::add);
            return lines;
        }

        return Collections.emptyList();
    }

    private static String geneCoverageHeader(final List<Integer> depthBuckets)
    {
        StringJoiner joiner = new StringJoiner(TSV_DELIM);

        joiner.add(FLD_GENE_NAME);
        joiner.add(FLD_CHROMOSOME);
        joiner.add(FLD_POS_START);
        joiner.add(FLD_POS_END);
        joiner.add(FLD_MV_LIKELIHOOD);

        for(int bucket = 0; bucket < depthBuckets.size() - 1; ++bucket)
        {
            int depth = depthBuckets.get(bucket);
            int depthNext = depthBuckets.get(bucket + 1);

            if(depthNext == depth + 1)
            {
                joiner.add(format("DR_%d", depth));
            }
            else
            {
                joiner.add(format("DR_%d_%d", depth, depthNext - 1));
            }
        }

        int maxDepth = depthBuckets.get(depthBuckets.size() - 1);
        joiner.add(format("DR_%d", maxDepth));

        return joiner.toString();
    }

    public static final String FLD_EXON_RANK = "ExonRank";
    public static final String FLD_MEDIAN_DEPTH = "MedianDepth";
    public static final String FLD_MEAN_DEPTH = "MeanDepth";
    public static final String FLD_PERC_ABOVE_DEPTH = "PercAboveDepth";

    public static String exonCoverageHeader(final int[] exonDepthThresholds)
    {
        StringJoiner sj = new StringJoiner(TSV_DELIM);

        sj.add(FLD_GENE_NAME);
        sj.add(FLD_CHROMOSOME);
        sj.add(FLD_POS_START);
        sj.add(FLD_POS_END);
        sj.add(FLD_EXON_RANK);
        sj.add(FLD_MEDIAN_DEPTH);
        sj.add(FLD_MEAN_DEPTH);

        for(int threshold : exonDepthThresholds)
        {
            sj.add(format("%s_%d", FLD_PERC_ABOVE_DEPTH, threshold));
        }

        return sj.toString();
    }

    private static String toString(final GeneDepth depth)
    {
        StringJoiner joiner = new StringJoiner(TSV_DELIM);
        joiner.add(depth.Gene);
        joiner.add(depth.Chromosome);
        joiner.add(String.valueOf(depth.PosStart));
        joiner.add(String.valueOf(depth.PosEnd));
        joiner.add(format("%.4f", depth.MissedVariantLikelihood));

        for(int i : depth.DepthCounts)
        {
            joiner.add(String.valueOf(i));
        }

        return joiner.toString();
    }
}
