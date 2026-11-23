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
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.IntStream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

public final class GeneDepthFile
{
    private static final String FLD_DEPTH_RANGE_PREFIX = "DR_";

    public static final String FLD_EXON_RANK = "ExonRank";
    public static final String FLD_MEDIAN_DEPTH = "MedianDepth";
    public static final String FLD_MEAN_DEPTH = "MeanDepth";

    public static String generateGeneCoverageFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + BAM_METRICS_FILE_ID + ".gene_coverage.tsv";
    }

    public static String generateExonCoverageFilename(final String basePath, final String sample)
    {
        return checkAddDirSeparator(basePath) + sample + BAM_METRICS_FILE_ID + ".exon_coverage.tsv";
    }

    public static void writeGeneCoverage(final String filename, final List<GeneDepth> depths) throws IOException
    {
        Files.write(new File(filename).toPath(), geneCoverageLines(depths));
    }

    public static void writeExonCoverage(final String filename, final List<ExonDepth> depths) throws IOException
    {
        Files.write(new File(filename).toPath(), exonCoverageLines(depths));
    }

    public static class GeneDepthComparator implements Comparator<GeneDepth>
    {
        public int compare(final GeneDepth first, final GeneDepth second)
        {
            return first.Gene.compareTo(second.Gene);
        }
    }

    public static List<GeneDepth> readGeneCoverage(final String filename) throws IOException
    {
        return readDepths(filename, false);
    }

    public static List<ExonDepth> readExonCoverage(final String filename) throws IOException
    {
        return readDepths(filename, true).stream().map(x -> (ExonDepth) x).toList();
    }

    private static List<GeneDepth> readDepths(final String filename, boolean forExons) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(filename).toPath());

        String header = lines.get(0);
        lines.remove(0);

        Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

        int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
        int geneIndex = fieldsIndexMap.get(FLD_GENE_NAME);
        int posStartIndex = fieldsIndexMap.get(FLD_POS_START);
        int posEndIndex = fieldsIndexMap.get(FLD_POS_END);

        int exonRankIndex = -1;
        int medianDepthIndex = -1;
        int meanDepthIndex = -1;

        if(forExons)
        {
            exonRankIndex = fieldsIndexMap.get(FLD_EXON_RANK);
            medianDepthIndex = fieldsIndexMap.get(FLD_MEDIAN_DEPTH);
            meanDepthIndex = fieldsIndexMap.get(FLD_MEAN_DEPTH);
        }

        String[] columns = header.split(TSV_DELIM, -1);
        int depthStartIndex = IntStream.range(0, columns.length)
                .filter(i -> columns[i].startsWith(FLD_DEPTH_RANGE_PREFIX))
                .findFirst()
                .orElse(-1);

        List<DepthRange> depthRanges = parseDepthRanges(columns, depthStartIndex);

        List<GeneDepth> depths = Lists.newArrayList();

        for(String line : lines)
        {
            String[] values = line.split(TSV_DELIM, -1);

            Map<DepthRange,Integer> depthCounts = Maps.newLinkedHashMap();

            for(int i = 0; i < depthRanges.size(); ++i)
            {
                int colIndex = depthStartIndex + i;
                depthCounts.put(depthRanges.get(i), Integer.parseInt(values[colIndex]));
            }

            String gene = values[geneIndex];
            String chromosome = values[chrIndex];
            int posStart = Integer.parseInt(values[posStartIndex]);
            int posEnd = Integer.parseInt(values[posEndIndex]);

            if(forExons)
            {
                depths.add(new ExonDepth(
                        gene, chromosome, posStart, posEnd,
                        Integer.parseInt(values[exonRankIndex]),
                        Double.parseDouble(values[medianDepthIndex]),
                        Double.parseDouble(values[meanDepthIndex]),
                        depthCounts));
            }
            else
            {
                depths.add(new GeneDepth(gene, chromosome, posStart, posEnd, depthCounts));
            }
        }

        return depths;
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

    private static List<String> geneCoverageLines(final List<GeneDepth> depths)
    {
        if(!depths.isEmpty())
        {
            depths.sort(new GeneDepthComparator());
            final List<String> lines = Lists.newArrayList();
            lines.add(geneCoverageHeader(depths.get(0).DepthCounts.keySet()));
            depths.stream().map(x -> toString(x, false)).forEach(lines::add);
            return lines;
        }

        return Collections.emptyList();
    }

    private static List<String> exonCoverageLines(final List<ExonDepth> depths)
    {
        if(depths.isEmpty())
            return Collections.emptyList();

        depths.sort(Comparator.comparing((ExonDepth x) -> x.Gene).thenComparingInt(x -> x.ExonRank));
        List<String> lines = Lists.newArrayList();
        lines.add(exonCoverageHeader(depths.get(0).DepthCounts.keySet()));
        depths.stream().map(x -> toString(x, true)).forEach(lines::add);
        return lines;
    }

    private static String geneCoverageHeader(final Collection<DepthRange> depthRanges)
    {
        return coverageHeader(depthRanges, false);
    }

    private static String exonCoverageHeader(final Collection<DepthRange> depthRanges)
    {
        return coverageHeader(depthRanges, true);
    }

    private static String coverageHeader(final Collection<DepthRange> depthRanges, boolean forExons)
    {
        StringJoiner joiner = new StringJoiner(TSV_DELIM);

        joiner.add(FLD_GENE_NAME);
        joiner.add(FLD_CHROMOSOME);
        joiner.add(FLD_POS_START);
        joiner.add(FLD_POS_END);

        if(forExons)
        {
            joiner.add(FLD_EXON_RANK);
            joiner.add(FLD_MEDIAN_DEPTH);
            joiner.add(FLD_MEAN_DEPTH);
        }

        for(DepthRange depthRange : depthRanges)
        {
            joiner.add(depthRangeFieldName(depthRange));
        }

        return joiner.toString();
    }

    public static String depthRangeFieldName(final DepthRange depthRange)
    {
        return depthRange.lower() == depthRange.upper()
                ? format("%s%d", FLD_DEPTH_RANGE_PREFIX, depthRange.lower())
                : format("%s%d_%d", FLD_DEPTH_RANGE_PREFIX, depthRange.lower(), depthRange.upper());
    }

    private static String toString(final GeneDepth depth, boolean forExons)
    {
        StringJoiner joiner = new StringJoiner(TSV_DELIM);
        joiner.add(depth.Gene);
        joiner.add(depth.Chromosome);
        joiner.add(String.valueOf(depth.PosStart));
        joiner.add(String.valueOf(depth.PosEnd));

        if(forExons)
        {
            ExonDepth exonDepth = (ExonDepth) depth;
            joiner.add(String.valueOf(exonDepth.ExonRank));
            joiner.add(format("%.0f", exonDepth.MedianDepth));
            joiner.add(format("%.2f", exonDepth.MeanDepth));
        }

        for(int count : depth.DepthCounts.values())
        {
            joiner.add(String.valueOf(count));
        }

        return joiner.toString();
    }
}
