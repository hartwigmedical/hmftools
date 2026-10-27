package com.hartwig.hmftools.common.metrics;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;

import org.junit.Test;

public class GeneDepthFileTest
{
    private static final Map<DepthRange,Integer> DEPTH_COUNTS = createDepthCounts();

    @Test
    public void canWriteAndReadGeneDepths() throws IOException
    {
        GeneDepth geneDepth = new GeneDepth("GENE_1", "1", 100, 200, DEPTH_COUNTS);

        File newFormatFile = createTempFile();
        GeneDepthFile.write(newFormatFile.getAbsolutePath(), Lists.newArrayList(geneDepth));

        List<String> expectedNewFormatLines = List.of(
                "GeneName\tChromosome\tPosStart\tPosEnd\tDR_29\tDR_30_39\tDR_100_149\tDR_10000",
                "GENE_1\t1\t100\t200\t10\t20\t30\t40"
        );
        assertEquals(expectedNewFormatLines, Files.readAllLines(newFormatFile.toPath()));

        List<String> oldFormatLines = List.of(
                "GeneName\tChromosome\tPosStart\tPosEnd\tMissedVariantLikelihood\tDR_29\tDR_30_39\tDR_100_149\tDR_10000",
                "GENE_1\t1\t100\t200\t0.0023\t10\t20\t30\t40"
        );
        File oldFormatFile = createTempFile();
        Files.write(oldFormatFile.toPath(), oldFormatLines);

        for(File file : List.of(newFormatFile, oldFormatFile))
        {
            List<GeneDepth> geneDepths = GeneDepthFile.read(file.getAbsolutePath());
            assertEquals(DEPTH_COUNTS, geneDepths.get(0).DepthCounts);
        }
    }

    private static Map<DepthRange,Integer> createDepthCounts()
    {
        Map<DepthRange,Integer> depthCounts = new LinkedHashMap<>();
        depthCounts.put(new DepthRange(29, 29), 10);
        depthCounts.put(new DepthRange(30, 39), 20);
        depthCounts.put(new DepthRange(100, 149), 30);
        depthCounts.put(new DepthRange(10000, 10000), 40);
        return depthCounts;
    }

    private static File createTempFile() throws IOException
    {
        File file = File.createTempFile("gene_coverage", ".tsv");
        file.deleteOnExit();
        return file;
    }
}
