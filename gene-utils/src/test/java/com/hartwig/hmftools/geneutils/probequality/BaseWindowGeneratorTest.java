package com.hartwig.hmftools.geneutils.probequality;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.junit.Test;

public class BaseWindowGeneratorTest
{
    private final BaseWindowGenerator mGenerator;

    public BaseWindowGeneratorTest() {
        MockRefGenome refGenome = new MockRefGenome();
        SpecificRegions specificRegions = new SpecificRegions();
        int baseWindowLength = 4;
        int baseWindowSpacing = 2;
        int batchSize = 3;
        mGenerator = new BaseWindowGenerator(refGenome, specificRegions, baseWindowLength, baseWindowSpacing, batchSize);
    }

    @Test
    public void test_batchRegions_single() {
        List<ChrBaseRegion> regions = List.of(
                new ChrBaseRegion("1", 1, 10),
                new ChrBaseRegion("2", 10, 20)
        );
        List<List<ChrBaseRegion>> expected = List.of(
                List.of(regions.get(0), regions.get(1))
        );
        List<List<ChrBaseRegion>> actual = mGenerator.batchRegions(regions.stream()).toList();
        assertEquals(expected, actual);
    }

    @Test
    public void test_batchRegions_multiple() {
        List<ChrBaseRegion> regions = List.of(
                new ChrBaseRegion("1", 1, 10),
                new ChrBaseRegion("2", 10, 20),
                new ChrBaseRegion("3", 20, 30),
                new ChrBaseRegion("4", 30, 40),
                new ChrBaseRegion("5", 40, 50),
                new ChrBaseRegion("6", 50, 60),
                new ChrBaseRegion("7", 60, 70)
        );
        List<List<ChrBaseRegion>> expected = List.of(
                List.of(regions.get(0), regions.get(1), regions.get(2)),
                List.of(regions.get(3), regions.get(4), regions.get(5)),
                List.of(regions.get(6))
        );
        List<List<ChrBaseRegion>> actual = mGenerator.batchRegions(regions.stream()).toList();
        assertEquals(expected, actual);
    }
}
