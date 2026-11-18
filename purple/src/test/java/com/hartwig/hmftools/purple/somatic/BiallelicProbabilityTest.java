package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.purple.somatic.BiallelicModel.calculateProbability;

import org.junit.Test;

import junit.framework.TestCase;

public class BiallelicProbabilityTest extends TestCase
{
    @Test
    public void testBiallelicCalculations()
    {
        double samplePurity = 0.34;
        int alleleReadCount = 869;
        int depth = 1932;
        double segmentCN = 4.04;
        double variantCN = 3.56;
        double minorAlleleCN = 2;

        double probability = calculateProbability(samplePurity, alleleReadCount, depth, segmentCN, variantCN, minorAlleleCN);

        assertEquals(0.7136, probability, 0.001);

        // scenario 1 - no LOH and not biallelic
        samplePurity = 1.0;
        alleleReadCount = 50;
        depth = 100;
        segmentCN = 2.0;
        variantCN = 1.0;
        minorAlleleCN = 1.0;

        probability = calculateProbability(samplePurity, alleleReadCount, depth, segmentCN, variantCN, minorAlleleCN);

        assertEquals(0, probability, 0.001);

        // scenario 1 - clear LOH and moderate biallelic likelihood
        samplePurity = 0.5;
        alleleReadCount = 25;
        depth = 100;
        segmentCN = 1.0;
        variantCN = 0.45;
        minorAlleleCN = 0.1;

        probability = calculateProbability(samplePurity, alleleReadCount, depth, segmentCN, variantCN, minorAlleleCN);

        assertEquals(0.373, probability, 0.001);

        /* old calc scenarios:

        // Scenario 3
        // input values:
        CN = 2.98;
        MACN = 0.991;
        VCN = 2.33;
        alleleReadCount = 36;

        // expected output:
        expectedBiallelicProbability = 0.054;

        correspondingBAF = 1 - (MACN / CN);
        copyNumber = createCopyNumber(CN, correspondingBAF);
        variant = createVariant(alleleReadCount, VCN);

        biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);

        // Scenario 4
        // input values:
        CN = 43.78;
        MACN = 0.0;
        VCN = 42.74;
        alleleReadCount = 627;

        // expected output:
        expectedBiallelicProbability = 0.442;

        correspondingBAF = 1 - (MACN / CN);
        copyNumber = createCopyNumber(CN, correspondingBAF);
        variant = createVariant(alleleReadCount, VCN);

        biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);

        // Scenario 5 - initially called BIALLELIC in COLO829v003T chr1:50084056
        // input values:
        CN = 2.81;
        MACN = 0.811;
        VCN = 2.83;
        alleleReadCount = 73;

        // expected output:
        expectedBiallelicProbability = 0.833;

        correspondingBAF = 1 - (MACN / CN);
        copyNumber = createCopyNumber(CN, correspondingBAF);
        variant = createVariant(alleleReadCount, VCN);

        biallelicProbability = SomaticPurityEnrichment.calculateBiallelic(copyNumber, variant);
        assertEquals(expectedBiallelicProbability, biallelicProbability, 0.01);


        */
    }
}