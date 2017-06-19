package com.hartwig.hmftools.common.purple.region;

import static com.hartwig.hmftools.common.purple.region.FittedRegionFactory.NORMAL_BAF;
import static com.hartwig.hmftools.common.purple.region.FittedRegionFactory.bafDeviation;
import static com.hartwig.hmftools.common.purple.region.FittedRegionFactory.cnvDeviation;
import static com.hartwig.hmftools.common.purple.region.FittedRegionFactory.modelBAF;
import static com.hartwig.hmftools.common.purple.region.FittedRegionFactory.modelBAFToMinimizeDeviation;
import static com.hartwig.hmftools.common.purple.region.FittedRegionFactory.modelRatio;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FittedRegionFactoryTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void expectedFit() {
        final FittedRegionFactory victim = new FittedRegionFactory(Gender.MALE, 12, 0.2);
        final FittedRegion result = victim.fitRegion(0.8, 0.7, create(180d / 280d + 0.01, 0.98 - 0.01));
        assertEquals(3, result.fittedPloidy());
        assertEquals(0.01, result.bafDeviation(), EPSILON);
        assertEquals(0.002, result.cnvDeviation(), EPSILON);
        assertEquals(0.011058009947947143, result.deviation(), EPSILON);
    }

    @NotNull
    private static ObservedRegion create(final double baf, final double ratio) {
        return ImmutableEnrichedRegion.builder()
                .observedBAF(baf)
                .bafCount(1)
                .chromosome("1")
                .start(1)
                .end(2)
                .observedTumorRatio(ratio)
                .observedNormalRatio(1)
                .build();
    }

    @Test
    public void diploidModelCNVRatio() {
        diploidModelCNVRatio(0.5, 0.4, 0.5);
        diploidModelCNVRatio(1, 0.4, 1);
        diploidModelCNVRatio(1, 0.6, 1);
    }

    @Test
    public void pureModelCNVRatio() {
        pureModelCNVRatio(1, 1, 2);
        pureModelCNVRatio(1.5, 1, 3);
        pureModelCNVRatio(2, 1, 4);

        pureModelCNVRatio(0.5, 0.5, 2);
        pureModelCNVRatio(0.75, 0.5, 3);
        pureModelCNVRatio(1, 0.5, 4);
    }

    @Test
    public void testModelCNVRatio() {
        assertModelCNVRatio(1.4, 0.8, 1, 3);
        assertModelCNVRatio(1.8, 0.8, 1, 4);

        assertModelCNVRatio(0.8 + 0.36, 0.9, 0.8, 3);
        assertModelCNVRatio(0.8 + 0.72, 0.9, 0.8, 4);

        assertModelCNVRatio(1.2, 1, 0.8, 3);
        assertModelCNVRatio(1.6, 1, 0.8, 4);
        assertModelCNVRatio(2.0, 1, 0.8, 5);
    }

    @Test
    public void testCNVDeviation() {
        assertCNVDeviation(0, 1, 0.3, 0.3);
        assertCNVDeviation(0, 0.8, 0.4, 0.4);

        assertCNVDeviation(0.05, 0.5, 1, 1.1);

        assertCNVDeviation(0.1, 1, 1, 1.1);
        assertCNVDeviation(0.1, 1, 1.5, 1.6);
        assertCNVDeviation(0.1, 1, 2, 2.1);
    }

    @Test
    public void testModelBAF() {
        assertModelBAF(18d / 28d, 0.8, 3, 2);
        assertModelBAF(26d / 28d, 0.8, 3, 3);

        assertModelBAF(NORMAL_BAF, 0.8, 4, 2);
    }

    @Test
    public void testModelBAFToMinimizeDeviation() {
        assertModelBAFToMinimizeDeviation(1, 1, 2, 1);

        assertModelBAFToMinimizeDeviation(18d / 28d, 0.8, 3, 0.65);
        assertModelBAFToMinimizeDeviation(26d / 28d, 0.8, 3, 0.95);

        assertModelBAFToMinimizeDeviation(17d / 27d, 0.7, 3, 0.65);
        assertModelBAFToMinimizeDeviation(24d / 27d, 0.7, 3, 0.95);
    }

    @Test
    public void testBAFDeviation() {
        assertEquals(0.01, bafDeviation(NORMAL_BAF, NORMAL_BAF + 0.01), EPSILON);
        assertEquals(0, bafDeviation(NORMAL_BAF, NORMAL_BAF), EPSILON);
        assertEquals(0.01, bafDeviation(NORMAL_BAF, NORMAL_BAF - 0.01), EPSILON);
        assertEquals(0.02, bafDeviation(NORMAL_BAF, NORMAL_BAF - 0.02), EPSILON);
    }

    private static void assertModelBAFToMinimizeDeviation(double expectedBAF, double purity, int ploidy,
            double actualBAF) {
        assertEquals(expectedBAF, modelBAFToMinimizeDeviation(purity, ploidy, actualBAF)[0], EPSILON);
    }

    private static void assertModelBAF(double expectedBAF, double purity, int ploidy, int betaAllele) {
        assertEquals(expectedBAF, modelBAF(purity, ploidy, betaAllele), EPSILON);
    }

    private static void assertCNVDeviation(double expectedDeviation, double cnvRatioWeighFactor, double modelCNVRatio,
            double tumorRatio) {
        assertEquals(expectedDeviation, cnvDeviation(cnvRatioWeighFactor, modelCNVRatio, tumorRatio), EPSILON);
    }

    private static void assertModelCNVRatio(double expectedRatio, double purity, double normFactor, int ploidy) {
        assertEquals(expectedRatio, modelRatio(purity, normFactor, ploidy), EPSILON);
    }

    private static void diploidModelCNVRatio(double expectedRatio, double purity, double normFactor) {
        assertEquals(expectedRatio, modelRatio(purity, normFactor, 2), EPSILON);
    }

    private static void pureModelCNVRatio(double expectedRatio, double normFactor, int ploidy) {
        assertEquals(expectedRatio, modelRatio(1, normFactor, ploidy), EPSILON);
    }
}
