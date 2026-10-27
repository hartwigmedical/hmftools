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
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public final class GeneDepthFile
{
    private static final String FLD_DEPTH_RANGE_PREFIX = "DR_";

    public static String generateGeneCoverageFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + BAM_METRICS_FILE_ID + ".gene_coverage.tsv";
    }

    public static String generateExonCoverageFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + BAM_METRICS_FILE_ID + ".exon_coverage.tsv";
    }

    public static void write(final String filename, final List<GeneDepth> depths) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(depths));
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
        String[] columns = header.split(TSV_DELIM, -1);
        int depthStartIndex = IntStream.range(0, columns.length)
                .filter(i -> columns[i].startsWith(FLD_DEPTH_RANGE_PREFIX))
                .findFirst()
                .orElse(-1);

        List<DepthRange> depthRanges = parseDepthRanges(columns, depthStartIndex);

        List<GeneDepth> geneDepths = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            Map<DepthRange,Integer> depthCounts = Maps.newLinkedHashMap();

            for(int i = 0; i < depthRanges.size(); ++i)
            {
                int colIndex = depthStartIndex + i;
                depthCounts.put(depthRanges.get(i), Integer.parseInt(values[colIndex]));
            }

            GeneDepth geneDepth = new GeneDepth(
                    values[geneIndex],
                    values[chrIndex],
                    Integer.parseInt(values[posStartIndex]),
                    Integer.parseInt(values[posEndIndex]),
                    depthCounts
            );

            geneDepths.add(geneDepth);
        }

        return geneDepths;
    }

    private static List<DepthRange> parseDepthRanges(final String[] columns, final int depthStartIndex)
    {
        List<DepthRange> depthRanges = Lists.newArrayList();

        for(int i = depthStartIndex; i < columns.length; ++i)
        {
            String depthRange = columns[i].substring(FLD_DEPTH_RANGE_PREFIX.length());
            String[] depthBounds = depthRange.split("_");

            int lowerDepth = Integer.parseInt(depthBounds[0]);
            int upperDepth = Integer.parseInt(depthBounds[depthBounds.length - 1]);

            depthRanges.add(new DepthRange(lowerDepth, upperDepth));
        }

        return depthRanges;
    }

    private static List<String> toLines(final List<GeneDepth> depths)
    {
        if(!depths.isEmpty())
        {
            Collections.sort(depths, new GeneDepthComparator());
            final List<String> lines = Lists.newArrayList();
            lines.add(geneCoverageHeader(depths.get(0).DepthCounts.keySet()));
            depths.stream().map(GeneDepthFile::toString).forEach(lines::add);
            return lines;
        }

        return Collections.emptyList();
    }

    private static String geneCoverageHeader(final Set<DepthRange> depthRanges)
    {
        StringJoiner joiner = new StringJoiner(TSV_DELIM);

        joiner.add(FLD_GENE_NAME);
        joiner.add(FLD_CHROMOSOME);
        joiner.add(FLD_POS_START);
        joiner.add(FLD_POS_END);

        for(DepthRange depthRange : depthRanges)
        {
            if(depthRange.lower() == depthRange.upper())
            {
                joiner.add(format("%s%d", FLD_DEPTH_RANGE_PREFIX, depthRange.lower()));
            }
            else
            {
                joiner.add(format("%s%d_%d", FLD_DEPTH_RANGE_PREFIX, depthRange.lower(), depthRange.upper()));
            }
        }

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

        for(int count : depth.DepthCounts.values())
        {
            joiner.add(String.valueOf(count));
        }

        return joiner.toString();
    }
}
