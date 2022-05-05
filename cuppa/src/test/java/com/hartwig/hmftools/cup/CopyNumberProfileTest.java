package com.hartwig.hmftools.cup;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.UNKNOWN;
import static com.hartwig.hmftools.cup.somatics.CopyNumberProfile.extractCopyNumberProfile;

import static junit.framework.TestCase.assertEquals;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.sigs.PositionFrequencies;

import org.junit.Test;

public class CopyNumberProfileTest
{
    @Test
    public void testBuildCopyNumberProfile()
    {
        List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();

        int bucketSize = 1000;
        int chrLength = 5000;
        Map<String,Integer> chrLengths = Maps.newHashMap();
        chrLengths.put("1", chrLength);
        chrLengths.put("2", chrLength);
        chrLengths.put("3", chrLength);
        PositionFrequencies posFrequencies = new PositionFrequencies(bucketSize, 1, chrLengths, false);

        assertEquals(15, posFrequencies.getBucketCount());

        String chromosome = "1";
        int cnRegionStart = 1;
        int cnRegionEnd = chrLength;
        copyNumbers.add(createCopyNumber(chromosome, cnRegionStart, cnRegionEnd, TELOMERE, 2));

        chromosome = "2";

        // initial region covers a few buckets
        cnRegionStart = 1;
        cnRegionEnd = (int)(bucketSize * 3.5);
        copyNumbers.add(createCopyNumber(chromosome, cnRegionStart, cnRegionEnd, UNKNOWN, 1));

        // then a few small ones within a bucket

        cnRegionStart = cnRegionEnd + 1;
        cnRegionEnd += 100;
        copyNumbers.add(createCopyNumber(chromosome, cnRegionStart, cnRegionEnd, UNKNOWN, 2));

        cnRegionStart = cnRegionEnd + 1;
        cnRegionEnd += 150;
        copyNumbers.add(createCopyNumber(chromosome, cnRegionStart, cnRegionEnd, UNKNOWN, 2));

        // then another straddling one
        cnRegionStart = cnRegionEnd + 1;
        cnRegionEnd = 4500;
        copyNumbers.add(createCopyNumber(chromosome, cnRegionStart, cnRegionEnd, UNKNOWN, 3));

        // then out to telomere
        cnRegionStart = cnRegionEnd + 1;
        copyNumbers.add(createCopyNumber(chromosome, cnRegionStart, chrLength, TELOMERE, 1));

        chromosome = "3";
        cnRegionStart = 1;
        cnRegionEnd = 250;
        copyNumbers.add(createCopyNumber(chromosome, cnRegionStart, cnRegionEnd, UNKNOWN, 2));

        cnRegionStart = cnRegionEnd + 1;
        cnRegionEnd = 750;
        copyNumbers.add(createCopyNumber(chromosome, cnRegionStart, cnRegionEnd, UNKNOWN, 1));

        cnRegionStart = cnRegionEnd + 1;
        cnRegionEnd = 2500;
        copyNumbers.add(createCopyNumber(chromosome, cnRegionStart, cnRegionEnd, UNKNOWN, 2));

        cnRegionStart = cnRegionEnd + 1;
        copyNumbers.add(createCopyNumber(chromosome, cnRegionStart, chrLength, TELOMERE, 4));

        double[] cnProfile = extractCopyNumberProfile(copyNumbers, posFrequencies);

        assertEquals(cnProfile.length, posFrequencies.getBucketCount());

        // chr 1
        assertEquals(cnProfile[0], 2.0, 0.1);
        assertEquals(cnProfile[4], 2.0, 0.1);

        // chr 2
        assertEquals(cnProfile[5], 1.0, 0.1);
        assertEquals(cnProfile[6], 1.0, 0.1);
        assertEquals(cnProfile[7], 1.0, 0.1);
        assertEquals(cnProfile[8], 1.75, 0.1);
        assertEquals(cnProfile[9], 2.75, 0.1);

        // chr 3
        assertEquals(cnProfile[10], 1.5, 0.1);
        assertEquals(cnProfile[11], 2.0, 0.1);
        assertEquals(cnProfile[12], 3.0, 0.1);
        assertEquals(cnProfile[13], 4.0, 0.1);
        assertEquals(cnProfile[14], 4.0, 0.1);
    }

    private static PurpleCopyNumber createCopyNumber(
            final String chromosome, int posStart, int posEnd, SegmentSupport endSupport, double cn)
    {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(chromosome)
                .start(posStart)
                .end(posEnd)
                .averageTumorCopyNumber(cn)
                .segmentStartSupport(SegmentSupport.NONE)
                .segmentEndSupport(endSupport)
                .method(CopyNumberMethod.UNKNOWN)
                .bafCount(0)
                .depthWindowCount(1)
                .gcContent(0)
                .minStart(0)
                .maxStart(0)
                .averageObservedBAF(0.5)
                .averageActualBAF(0.5)
                .build();
    }


}
