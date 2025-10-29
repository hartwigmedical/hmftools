package com.hartwig.hmftools.common.metrics;

import static com.hartwig.hmftools.common.metrics.BamMetricsCommon.BAM_METRICS_FILE_ID;
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
    public static final String COL_GENE = "gene";
    public static final String COL_CHROMOSOME = "chromosome";
    public static final String COL_POS_START = "posStart";
    public static final String COL_POS_END = "posEnd";
    public static final String COL_MV_LIKELIHOOD = "missedVariantLikelihood";

    public static String generateGeneCoverageFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + BAM_METRICS_FILE_ID + ".gene_coverage.tsv";
    }

    public static String generateExonMediansFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + BAM_METRICS_FILE_ID + ".exon_medians.tsv";
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

    public static List<int[]> readDepthRanges(final String filename) throws IOException
    {
        List<int[]> ranges = Lists.newArrayList();

        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String header = lines.get(0);
        String[] values = header.split(TSV_DELIM, -1);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);
        int mvlIndex = fieldsIndexMap.get(COL_MV_LIKELIHOOD);
        int depthStartIndex = fieldsIndexMap.get(COL_MV_LIKELIHOOD) + 1;
        int depthColumns = values.length - depthStartIndex;

        for(int i = 0; i < depthColumns; ++i)
        {
            int colIndex = depthStartIndex + i;
            String colStr = values[colIndex];

            // depth columns are like: DR_0 or DR_4000_4999
            String[] depthParts = colStr.split("_");
            int rangeStart = Integer.parseInt(depthParts[1]);
            int rangeEnd = depthParts.length == 3 ? Integer.parseInt(depthParts[2]) : rangeStart;
            ranges.add(new int[] {rangeStart, rangeEnd});
        }

        return ranges;
    }

    public static List<GeneDepth> read(final String filename) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String header = lines.get(0);
        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        int chrIndex = fieldsIndexMap.get(COL_CHROMOSOME);
        int geneIndex = fieldsIndexMap.get(COL_GENE);
        int posStartIndex = fieldsIndexMap.get(COL_POS_START);
        int posEndIndex = fieldsIndexMap.get(COL_POS_END);
        int mvlIndex = fieldsIndexMap.get(COL_MV_LIKELIHOOD);
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
            lines.add(header(depthBuckets));
            depths.stream().map(GeneDepthFile::toString).forEach(lines::add);
            return lines;
        }

        return Collections.emptyList();
    }

    private static String header(final List<Integer> depthBuckets)
    {
        StringJoiner joiner = new StringJoiner(TSV_DELIM);

        joiner.add(COL_GENE);
        joiner.add(COL_CHROMOSOME);
        joiner.add(COL_POS_START);
        joiner.add(COL_POS_END);
        joiner.add(COL_MV_LIKELIHOOD);

        for(int bucket = 0; bucket < depthBuckets.size() - 1; ++bucket)
        {
            int depth = depthBuckets.get(bucket);
            int depthNext = depthBuckets.get(bucket + 1);

            if(depthNext == depth + 1)
            {
                joiner.add(String.format("DR_%d", depth));
            }
            else
            {
                joiner.add(String.format("DR_%d_%d", depth, depthNext - 1));
            }
        }

        int maxDepth = depthBuckets.get(depthBuckets.size() - 1);
        joiner.add(String.format("DR_%d", maxDepth));

        return joiner.toString();
    }

    private static String toString(final GeneDepth depth)
    {
        StringJoiner joiner = new StringJoiner(TSV_DELIM);
        joiner.add(depth.Gene);
        joiner.add(depth.Chromosome);
        joiner.add(String.valueOf(depth.PosStart));
        joiner.add(String.valueOf(depth.PosEnd));
        joiner.add(String.format("%.4f", depth.MissedVariantLikelihood));

        for(int i : depth.DepthCounts)
        {
            joiner.add(String.valueOf(i));
        }

        return joiner.toString();
    }
}
