package com.hartwig.hmftools.common.purple.copynumber;

import static org.apache.commons.math3.util.Precision.EPSILON;
import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.sv.ImmutableStructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class StructuralVariantImpliedTest {

    private static final String CHROMOSOME = "1";
    private static final int PLOIDY = 1;
    private static final PurityAdjuster PURE = new PurityAdjuster(Gender.FEMALE, 1, 1);

    @Test
    public void testNonSymmetricMultiPass() {
        final StructuralVariant firstSV = sv(1001, 4001, StructuralVariantType.DEL, 0.25, 0.25);
        final StructuralVariant secondSV = sv(2001, 3001, StructuralVariantType.DEL, 1 / 3d, 1 / 3d);

        final CombinedRegion firstCN = copyNumber(1, 1000, 40, SegmentSupport.NONE);
        final CombinedRegion secondCN = copyNumber(1001, 2000, 0, SegmentSupport.DEL);
        final CombinedRegion thirdCN = copyNumber(2001, 3000, 0, SegmentSupport.DEL);
        final CombinedRegion forthCN = copyNumber(3001, 4000, 0, SegmentSupport.DEL);
        final CombinedRegion fifthCN = copyNumber(4001, 5000, 10, SegmentSupport.NONE);

        final List<StructuralVariant> svs = Lists.newArrayList(firstSV, secondSV);
        final ListMultimap<String, CombinedRegion> copyNumbers = copyNumbers(firstCN, secondCN, thirdCN, forthCN, fifthCN);

        final StructuralVariantImplied victim = new StructuralVariantImplied(PURE);
        final List<CombinedRegion> result = victim.svImpliedCopyNumber(svs, copyNumbers).get(CHROMOSOME);
        assertEquals(5, result.size());
        assertEquals(40.00, result.get(0).tumorCopyNumber(), EPSILON);
        assertEquals(33.75, result.get(1).tumorCopyNumber(), EPSILON);
        assertEquals(12.50, result.get(2).tumorCopyNumber(), EPSILON);
        assertEquals(03.75, result.get(3).tumorCopyNumber(), EPSILON);
        assertEquals(10.00, result.get(4).tumorCopyNumber(), EPSILON);
    }

    @Test
    public void testSVSmoothing() {
        final StructuralVariant firstSV = sv(1001, 4001, StructuralVariantType.DEL, 0.25, 0.25);
        final StructuralVariant secondSV = sv(2001, 3001, StructuralVariantType.DEL, 1 / 3d, 1 / 3d);

        final CombinedRegion firstCN = copyNumber(1, 1000, 40, SegmentSupport.NONE);
        final CombinedRegion secondCN = copyNumber(1001, 2000, 0, SegmentSupport.DEL);
        final CombinedRegion thirdCN = copyNumber(2001, 3000, 0, SegmentSupport.NONE);
        final CombinedRegion forthCN = copyNumber(3001, 4000, 0, SegmentSupport.DEL);
        final CombinedRegion fifthCN = copyNumber(4001, 5000, 10, SegmentSupport.NONE);

        final List<StructuralVariant> svs = Lists.newArrayList(firstSV, secondSV);
        final ListMultimap<String, CombinedRegion> copyNumbers = copyNumbers(firstCN, secondCN, thirdCN, forthCN, fifthCN);

        final StructuralVariantImplied victim = new StructuralVariantImplied(PURE);
        final List<CombinedRegion> result = victim.svImpliedCopyNumber(svs, copyNumbers).get(CHROMOSOME);
        assertEquals(4, result.size());
        assertEquals(40.00, result.get(0).tumorCopyNumber(), EPSILON);
        assertEquals(33.75, result.get(1).tumorCopyNumber(), EPSILON);
        assertEquals(03.75, result.get(2).tumorCopyNumber(), EPSILON);
        assertEquals(10.00, result.get(3).tumorCopyNumber(), EPSILON);
    }

    @Test
    public void testImpliedCopyNumber() {
        final StructuralVariantLegPloidy left = create(1, Optional.of(4d), Optional.empty());
        final StructuralVariantLegPloidy right = create(-1, Optional.empty(), Optional.of(5d));

        final double bothKnown = StructuralVariantImplied.inferCopyNumberFromStructuralVariants(Optional.of(left), Optional.of(right));
        assertEquals(3.5, bothKnown, EPSILON);

        final double leftKnown = StructuralVariantImplied.inferCopyNumberFromStructuralVariants(Optional.of(left), Optional.empty());
        assertEquals(3, leftKnown, EPSILON);

        final double rightKnown = StructuralVariantImplied.inferCopyNumberFromStructuralVariants(Optional.empty(), Optional.of(right));
        assertEquals(4, rightKnown, EPSILON);
    }

    @NotNull
    private static StructuralVariant sv(long start, long end, StructuralVariantType type, double startAF, double endAF) {
        return PurpleDatamodelTest.createStructuralVariant(CHROMOSOME, start, CHROMOSOME, end, type, startAF, endAF).build();
    }

    @NotNull
    private static CombinedRegion copyNumber(long start, long end, double copyNumber, SegmentSupport support) {
        final FittedRegion region = PurpleDatamodelTest.createDefaultFittedRegion(CHROMOSOME, start, end)
                .tumorCopyNumber(copyNumber)
                .tumorBAF(0.5)
                .support(support)
                .build();

        final CombinedRegion result = new CombinedRegion(true, region);
        if (Doubles.positive(copyNumber)) {
            result.setCopyNumberMethod(CopyNumberMethod.BAF_WEIGHTED);
        }

        return result;
    }

    @NotNull
    private static ListMultimap<String, CombinedRegion> copyNumbers(CombinedRegion... copyNumbers) {
        final ListMultimap<String, CombinedRegion> result = ArrayListMultimap.create();
        result.putAll(CHROMOSOME, Lists.newArrayList(copyNumbers));
        return result;
    }

    @NotNull
    private static StructuralVariantLegPloidy create(int orientation, @NotNull final Optional<Double> leftCopyNumber,
            @NotNull final Optional<Double> rightCopyNumber) {
        return ImmutableStructuralVariantLegPloidy.builder()
                .chromosome(CHROMOSOME)
                .position(1)
                .orientation((byte) orientation)
                .vaf(0.5)
                .alleleFrequency(0.5)
                .homology("")
                .weight(1)
                .averageImpliedPloidy(PLOIDY)
                .unweightedImpliedPloidy(PLOIDY)
                .leftCopyNumber(leftCopyNumber)
                .rightCopyNumber(rightCopyNumber)
                .build();
    }
}
