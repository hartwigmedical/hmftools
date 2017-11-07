package com.hartwig.hmftools.purple.plot;

import static com.google.common.collect.Lists.newArrayList;
import static com.hartwig.hmftools.purple.plot.CopyNumberCharts.minorAllele;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.copynumber.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;
import org.jfree.data.xy.CategoryTableXYDataset;
import org.junit.Test;

public class CopyNumberChartsTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testNegativeMinorBafs() {
        CategoryTableXYDataset dataset = minorAllele(newArrayList(createCopyNumber(1.1, 1)));
        assertEquals(-0.1, dataset.getX(0, 0).doubleValue(), EPSILON);
        assertEquals(10, dataset.getY(0, 0).doubleValue(), EPSILON);

        dataset = minorAllele(newArrayList(createCopyNumber(1.7, 1)));
        assertEquals(-0.7, dataset.getX(0, 0).doubleValue(), EPSILON);
        assertEquals(10, dataset.getY(0, 0).doubleValue(), EPSILON);
    }

    @Test
    public void testLowerBoundary() {
        CategoryTableXYDataset dataset = minorAllele(newArrayList(createCopyNumber(2, 1)));
        assertEquals(-1, dataset.getX(0, 0).doubleValue(), EPSILON);
        assertEquals(10, dataset.getY(0, 0).doubleValue(), EPSILON);

        dataset = minorAllele(newArrayList(createCopyNumber(3, 1)));
        assertEquals(-1, dataset.getX(0, 0).doubleValue(), EPSILON);
        assertEquals(10, dataset.getY(0, 0).doubleValue(), EPSILON);
    }

    @Test
    public void testUpperBoundary() {
        CategoryTableXYDataset dataset = minorAllele(newArrayList(createCopyNumber(0, 5)));
        assertEquals(5, dataset.getX(0, 0).doubleValue(), EPSILON);
        assertEquals(10, dataset.getY(0, 0).doubleValue(), EPSILON);

        dataset = minorAllele(newArrayList(createCopyNumber(0, 6)));
        assertEquals(5, dataset.getX(0, 0).doubleValue(), EPSILON);
        assertEquals(10, dataset.getY(0, 0).doubleValue(), EPSILON);
    }

    @NotNull
    private static PurpleCopyNumber createCopyNumber(double baf, double copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome("1")
                .start(1)
                .end(100)
                .averageTumorCopyNumber(copyNumber)
                .bafCount(10)
                .averageObservedBAF(baf)
                .ratioSupport(true)
                .support(SegmentSupport.NONE)
                .averageActualBAF(baf).build();
    }
}
