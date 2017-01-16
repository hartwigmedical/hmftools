package com.hartwig.hmftools.patientreporter.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.Collection;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.patientreporter.slicing.GenomeRegion;

import org.junit.Test;

public class CopyNumberAnalyserTest {

    private static final String CHROMOSOME = "X";
    private static final double EPSILON = 1.0e-10;

    @Test
    public void worksAsExpected() {
        final GenomeRegion first = new GenomeRegion(CHROMOSOME, 101, 200);
        final GenomeRegion second = new GenomeRegion(CHROMOSOME, 301, 400);
        final GenomeRegion third = new GenomeRegion(CHROMOSOME, 601, 700);
        final Collection<GenomeRegion> regions = Lists.newArrayList(first, second, third);

        final Collection<CopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(new CopyNumber(CHROMOSOME, 11, 20, 1));
        copyNumbers.add(new CopyNumber(CHROMOSOME, 121, 170, 4));
        copyNumbers.add(new CopyNumber(CHROMOSOME, 191, 260, 5));
        copyNumbers.add(new CopyNumber(CHROMOSOME, 291, 500, 1));

        final Map<GenomeRegion, CopyNumberStats> stats = CopyNumberAnalyser.run(regions, copyNumbers);

        final CopyNumberStats firstStat = stats.get(first);
        assertEquals(2, firstStat.min());
        assertEquals(5, firstStat.max());
        assertEquals(3.3, firstStat.mean(), EPSILON);

        final CopyNumberStats secondStat = stats.get(second);
        assertEquals(1, secondStat.min());
        assertEquals(1, secondStat.max());
        assertEquals(1D, secondStat.mean(), EPSILON);

        final CopyNumberStats thirdStat = stats.get(third);
        assertEquals(2, thirdStat.min());
        assertEquals(2, thirdStat.max());
        assertEquals(2D, thirdStat.mean(), EPSILON);
    }
}