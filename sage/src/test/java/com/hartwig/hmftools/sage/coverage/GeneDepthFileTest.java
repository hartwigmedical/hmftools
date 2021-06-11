package com.hartwig.hmftools.sage.coverage;

import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class GeneDepthFileTest
{

    @Test
    public void testFormat()
    {
        int[] baseCoverage = new int[37];
        baseCoverage[15] = 100;
        baseCoverage[20] = 100;
        baseCoverage[31] = 800;

        GeneDepth geneDepth = ImmutableGeneDepth.builder()
                .gene("TP53")
                .depthCounts(baseCoverage)
                .missedVariantLikelihood(GeneCoverage.missedVariantLikelihood(baseCoverage))
                .build();

        String geneDepthString = GeneDepthFile.toString(geneDepth);
        assertTrue(geneDepthString.startsWith("TP53\t0.23%\t0"));
    }

}
