package com.hartwig.hmftools.common.sage;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

public final class GeneDepthFile
{
    public static final String DELIM = "\t";

    public static final String COL_GENE = "gene";
    public static final String COL_CHROMOSOME = "chromosome";
    public static final String COL_POS_START = "posStart";
    public static final String COL_POS_END = "posEnd";
    public static final String COL_MV_LIKELIHOOD = "missedVariantLikelihood";

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
        StringJoiner joiner = new StringJoiner(DELIM);

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
        StringJoiner joiner = new StringJoiner(DELIM);
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
