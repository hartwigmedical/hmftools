package com.hartwig.hmftools.purple.region;

import static com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod.BAF_WEIGHTED;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.BND;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.purple.tools.HrdDetection;

import org.junit.Test;

public class HrdDetectionTest
{
    private static final String CHR_1 = "1";
    private static final String CHR_2 = "2";

    @Test
    public void testLohSegments()
    {
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();

        int minLohLength = 20000;

        copyNumbers.add(makeCopyNumber(TELOMERE, BND, 1, 10000, 2, 0.5));

        // start LOH
        copyNumbers.add(makeCopyNumber(BND, BND, 10001, 20000, 1, 1.0));

        // continues
        copyNumbers.add(makeCopyNumber(BND, BND, 20001, 30000, 1, 1.0));

        // broken by short insert
        copyNumbers.add(makeCopyNumber(BND, BND, 30001, 30100, 2, 0.5));

        // continues
        copyNumbers.add(makeCopyNumber(BND, BND, 30101, 50000, 1, 1.0));

        // ends
        copyNumbers.add(makeCopyNumber(BND, BND, 50001, 60000, 4, 0.75));

        HrdDetection hrdDetection = new HrdDetection(10000, 10000, minLohLength, 5);

        assertEquals(1, hrdDetection.calcLohSegments(copyNumbers));

        // doesn't call short LOHs
        copyNumbers.clear();

        copyNumbers.add(makeCopyNumber(TELOMERE, BND, 1, 10000, 2, 0.5));
        copyNumbers.add(makeCopyNumber(BND, BND, 10001, 20000, 1, 1.0));
        copyNumbers.add(makeCopyNumber(BND, BND, 20001, 60000, 2, 0.5));

        assertEquals(0, hrdDetection.calcLohSegments(copyNumbers));

        // can call from and to telomere
        copyNumbers.clear();

        copyNumbers.add(makeCopyNumber(CHR_1, TELOMERE, CENTROMERE, 1, 100000, 1, 1));
        copyNumbers.add(makeCopyNumber(CHR_1, CENTROMERE, TELOMERE, 100001, 500000, 1, 1));

        // each side of the centromere
        copyNumbers.add(makeCopyNumber(CHR_2, TELOMERE, BND, 1, 10000, 4, 0.5));
        copyNumbers.add(makeCopyNumber(CHR_2, BND, CENTROMERE, 10001, 100000, 1, 1));
        copyNumbers.add(makeCopyNumber(CHR_2, CENTROMERE, BND, 100001, 200000, 1, 1));
        copyNumbers.add(makeCopyNumber(CHR_2, BND, TELOMERE, 200001, 500000, 4, 0.5));

        assertEquals(4, hrdDetection.calcLohSegments(copyNumbers));
    }

    @Test
    public void testImbalances()
    {
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();

        int imbalanceLength = 20000;

        // imbalance from the telomere
        copyNumbers.add(makeCopyNumber(CHR_1, TELOMERE, BND, 1, 30000, 3, 0.67));

        copyNumbers.add(makeCopyNumber(CHR_1, BND, CENTROMERE, 30001, 100000, 2, 0.5));
        copyNumbers.add(makeCopyNumber(CHR_1, CENTROMERE, BND, 100001, 110000, 2, 0.5));

        // another out to the telomere
        copyNumbers.add(makeCopyNumber(CHR_1, BND, TELOMERE, 110001, 200000, 1.75, 0.8));

        // same again but too short
        copyNumbers.add(makeCopyNumber(CHR_2, TELOMERE, BND, 1, 10000, 3, 0.67));

        copyNumbers.add(makeCopyNumber(CHR_2, BND, CENTROMERE, 10001, 100000, 2, 0.5));
        copyNumbers.add(makeCopyNumber(CHR_2, CENTROMERE, BND, 100001, 190000, 2, 0.5));

        // another out to the telomere
        copyNumbers.add(makeCopyNumber(CHR_2, BND, TELOMERE, 190001, 200000, 1.75, 0.8));

        HrdDetection hrdDetection = new HrdDetection(10000, imbalanceLength, 10000, 5);

        assertEquals(2, hrdDetection.calcSegmentImbalances(copyNumbers));

        copyNumbers.clear();

        // cannot run into the centromere
        copyNumbers.add(makeCopyNumber(CHR_1, TELOMERE, BND, 1, 90000, 1, 1));
        copyNumbers.add(makeCopyNumber(CHR_1, BND, BND, 90001, 95000, 1, 1));
        copyNumbers.add(makeCopyNumber(CHR_1, BND, CENTROMERE, 95001, 100000, 1, 1));
        copyNumbers.add(makeCopyNumber(CHR_1, CENTROMERE, TELOMERE, 100001, 200000, 1, 1));

        assertEquals(0, hrdDetection.calcSegmentImbalances(copyNumbers));

        // check ignores short TIs
        copyNumbers.clear();

        copyNumbers.add(makeCopyNumber(CHR_1, TELOMERE, CENTROMERE, 1, 100000, 2, 0.5));
        copyNumbers.add(makeCopyNumber(CHR_1, CENTROMERE, BND, 100001, 110000, 2, 0.5));
        copyNumbers.add(makeCopyNumber(CHR_1, BND, BND, 110001, 120000, 1, 1));
        copyNumbers.add(makeCopyNumber(CHR_1, BND, BND, 120001, 130000, 2, 0.5)); // stops
        copyNumbers.add(makeCopyNumber(CHR_1, BND, BND, 130001, 190000, 1, 1)); // starts again
        copyNumbers.add(makeCopyNumber(CHR_1, BND, BND, 190001, 190500, 2, 0.5)); // short TI doesn't interupt
        copyNumbers.add(makeCopyNumber(CHR_1, BND, TELOMERE, 190501, 200000, 1, 1)); // starts again

        assertEquals(1, hrdDetection.calcSegmentImbalances(copyNumbers));
    }

        private static PurpleCopyNumber makeCopyNumber(
            final SegmentSupport segStart, final SegmentSupport segEnd, int posStart, int posEnd,
            double copyNumber, double baf)
    {
        return makeCopyNumber(CHR_1, segStart, segEnd, posStart, posEnd, copyNumber, baf);
    }

    private static PurpleCopyNumber makeCopyNumber(
            final String chromosome, final SegmentSupport segStart, final SegmentSupport segEnd, int posStart, int posEnd,
            double copyNumber, double baf)
    {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(posStart)
                .end(posEnd)
                .bafCount(1)
                .averageTumorCopyNumber(copyNumber)
                .averageObservedBAF(baf)
                .averageActualBAF(baf)
                .segmentStartSupport(segStart)
                .segmentEndSupport(segEnd)
                .method(BAF_WEIGHTED)
                .gcContent(0.5)
                .depthWindowCount(1)
                .minStart(posStart)
                .maxStart(posEnd).build();
    }
}
