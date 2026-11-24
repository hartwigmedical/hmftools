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

        assertEquals(0.7333, probability, 0.001);

        // test 1 - extreme subclonality - relative calc uses fallback value giving zero biallelic prob
        samplePurity = 0.32;
        alleleReadCount = 8;
        depth = 755;
        segmentCN = 1.58;
        variantCN = 0.062;
        minorAlleleCN = 0.552;

        probability = calculateProbability(samplePurity, alleleReadCount, depth, segmentCN, variantCN, minorAlleleCN);

        assertEquals(0, probability, 0.001);

        // test 2 - depth dampening plus AF uncertainty drags lowers p(noWT|LOH) certainty, adjustment further drags down p(noWT|!LOH) estimate
        samplePurity = 0.73;
        alleleReadCount = 1337;
        depth = 2144;
        segmentCN = 2.74;
        variantCN = 2.17;
        minorAlleleCN = 0.889;

        probability = calculateProbability(samplePurity, alleleReadCount, depth, segmentCN, variantCN, minorAlleleCN);

        assertEquals(0.2989, probability, 0.001);

        // test 3 - adjustment to p(NOWT|!LOH) adjustment (since VCN is closer to MACN than CN) drags p(noWT|!LOH) down further
        samplePurity = 0.26;
        alleleReadCount = 11;
        depth = 107;
        segmentCN = 1.61;
        variantCN = 0.751;
        minorAlleleCN = 0.489;

        probability = calculateProbability(samplePurity, alleleReadCount, depth, segmentCN, variantCN, minorAlleleCN);

        assertEquals(0.0399, probability, 0.001);

        // test 4 - heavy adjustment to p(NOWT|!LOH) (because VCN is very close to CN and not close to either allele CN) drags the
        // biallelic prob way up
        samplePurity = 0.45;
        alleleReadCount = 20;
        depth = 58;
        segmentCN = 1.39;
        variantCN = 1.32;
        minorAlleleCN = 0.63;

        probability = calculateProbability(samplePurity, alleleReadCount, depth, segmentCN, variantCN, minorAlleleCN);

        assertEquals(0.9138, probability, 0.001);
    }
}