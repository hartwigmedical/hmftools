package com.hartwig.hmftools.sage.coverage;

import static com.hartwig.hmftools.sage.SageCommon.DELIM;
import static com.hartwig.hmftools.sage.coverage.GeneCoverage.DEPTH_BUCKETS;
import static com.hartwig.hmftools.sage.coverage.GeneCoverage.MAX_DEPTH_BUCKET;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class GeneDepthFile
{
    public static void write(@NotNull final String filename, @NotNull final List<GeneDepth> depths) throws IOException
    {
        Files.write(new File(filename).toPath(), toLines(depths));
    }

    @NotNull
    private static List<String> toLines(@NotNull final List<GeneDepth> depths)
    {
        if(!depths.isEmpty())
        {
            final List<String> lines = Lists.newArrayList();
            lines.add(header());
            depths.stream().map(GeneDepthFile::toString).forEach(lines::add);
            return lines;
        }

        return Collections.emptyList();
    }

    static String header()
    {
        StringJoiner joiner = new StringJoiner(DELIM);

        joiner.add("gene").add("missedVariantLikelihood");

        for(int bucket = 0; bucket < DEPTH_BUCKETS.size() - 1; ++bucket)
        {
            int depth = DEPTH_BUCKETS.get(bucket);
            int depthNext = DEPTH_BUCKETS.get(bucket + 1);

            if(depthNext == depth + 1)
            {
                joiner.add(String.format("DR_%d", depth));
            }
            else
            {
                joiner.add(String.format("DR_%d_%d", depth, depthNext - 1));
            }
        }

        joiner.add(String.format("DR_%d", MAX_DEPTH_BUCKET));

        return joiner.toString();
    }

    @NotNull
    static String toString(@NotNull final GeneDepth depth)
    {
        StringJoiner joiner = new StringJoiner(DELIM);
        joiner.add(depth.Gene);
        joiner.add(String.format("%.4f", depth.MissedVariantLikelihood));

        for(int i : depth.DepthCounts)
        {
            joiner.add(String.valueOf(i));
        }

        return joiner.toString();
    }

}
