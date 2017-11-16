package com.hartwig.hmftools.common.purple.copynumber.sv;

import static java.util.Collections.singleton;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotEquals;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class StructuralVariantPloidyFactoryTest {

    private static final String CHROMOSOME = "X";
    private static final double EPSILON = 1e-10;
    private static final PurityAdjuster PURE = new PurityAdjuster(Gender.FEMALE, 1d, 1d);
    private static final StructuralVariantPloidyFactory PURE_PLOIDY_FACTORY = new StructuralVariantPloidyFactory(PURE);

    @Test
    public void testSingleValidLeg() {
        final StructuralVariantLeg start = leg(1001, 1, 0.25);

        final PurpleCopyNumber left = copyNumber(1, 1000, 4);
        final PurpleCopyNumber middle = copyNumber(1001, 2000, 3);
        final PurpleCopyNumber right = copyNumber(2001, 3000, 4);

        final StructuralVariantLegs legs = ImmutableStructuralVariantLegs.builder().start(start).end(Optional.empty()).build();
        final ListMultimap<String, PurpleCopyNumber> copyNumbers = copyNumbers(left, middle, right);
        final List<StructuralVariantPloidy> ploidies = PURE_PLOIDY_FACTORY.create(legs, copyNumbers);
        assertEquals(1, ploidies.size());
        for (StructuralVariantPloidy ploidy : ploidies) {
            assertEquals(1d, ploidy.averageImpliedPloidy(), EPSILON);
            assertEquals(1d, ploidy.weight(), EPSILON);
        }
    }

    @Test
    public void testBothLegsValid() {
        final StructuralVariantLeg start = leg(1001, 1, 0.25);
        final StructuralVariantLeg end = leg(2001, -1, 0.25);

        final PurpleCopyNumber left = copyNumber(1, 1000, 4);
        final PurpleCopyNumber middle = copyNumber(1001, 2000, 3);
        final PurpleCopyNumber right = copyNumber(2001, 3000, 4);

        final StructuralVariantLegs legs = ImmutableStructuralVariantLegs.builder().start(start).end(end).build();
        final ListMultimap<String, PurpleCopyNumber> copyNumbers = copyNumbers(left, middle, right);
        final List<StructuralVariantPloidy> ploidies = PURE_PLOIDY_FACTORY.create(legs, copyNumbers);
        assertEquals(2, ploidies.size());
        for (StructuralVariantPloidy ploidy : ploidies) {
            assertEquals(1d, ploidy.averageImpliedPloidy(), EPSILON);
            assertEquals(2d, ploidy.weight(), EPSILON);
        }
    }

    @Test
    public void testSelectCorrectOrAlternativeCopyNumberForLeg() {
        final StructuralVariantPloidyFactory pureFactory = new StructuralVariantPloidyFactory(PURE);

        final StructuralVariantLeg positiveLeg = leg(1001, 1, 0.25);
        final StructuralVariantLeg negativeLeg = leg(2001, -1, 0.25);
        final PurpleCopyNumber left = copyNumber(1, 1000, 4);
        final PurpleCopyNumber middle = copyNumber(1001, 2000, 3);
        final PurpleCopyNumber right = copyNumber(2001, 3000, 4);

        assertPloidy(1, false, pureFactory.create(positiveLeg, GenomeRegionSelectorFactory.create(singleton(left))));
        assertPloidy(1, true, pureFactory.create(positiveLeg, GenomeRegionSelectorFactory.create(singleton(middle))));

        assertPloidy(1, true, pureFactory.create(negativeLeg, GenomeRegionSelectorFactory.create(singleton(middle))));
        assertPloidy(1, false, pureFactory.create(negativeLeg, GenomeRegionSelectorFactory.create(singleton(right))));
    }

    @Test
    public void testPurityAdjustedPloidy() {

        final StructuralVariantLeg leg = leg(1001, 1, 0.5);
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList(copyNumber(1, 1000, 2), copyNumber(1001, 200, 1));

        final StructuralVariantPloidyFactory pureFactory = new StructuralVariantPloidyFactory(PURE);
        Optional<ModifiableStructuralVariantPloidy> purePloidy = pureFactory.create(leg, GenomeRegionSelectorFactory.create(copyNumbers));
        assertPloidy(1d, purePloidy);

        final PurityAdjuster diluted = new PurityAdjuster(Gender.FEMALE, 0.8, 1);
        final StructuralVariantPloidyFactory dilutedFactory = new StructuralVariantPloidyFactory(diluted);
        Optional<ModifiableStructuralVariantPloidy> dilutedPloidy =
                dilutedFactory.create(leg, GenomeRegionSelectorFactory.create(copyNumbers));
        assertPloidy(1.25d, dilutedPloidy);

        final PurityAdjuster male = new PurityAdjuster(Gender.MALE, 0.8, 1);
        final StructuralVariantPloidyFactory maleFactory = new StructuralVariantPloidyFactory(male);
        Optional<ModifiableStructuralVariantPloidy> malePloidy = maleFactory.create(leg, GenomeRegionSelectorFactory.create(copyNumbers));
        assertPloidy(1.125d, malePloidy);
    }

    private static StructuralVariantLeg leg(long position, int orientation, double vaf) {
        return ImmutableStructuralVariantLeg.builder().chromosome(CHROMOSOME).position(position).orientation(orientation).vaf(vaf).build();
    }

    private void assertPloidy(double expected, @NotNull final Optional<ModifiableStructuralVariantPloidy> ploidy) {
        assertEquals(expected, ploidy.map(ModifiableStructuralVariantPloidy::unweightedImpliedPloidy).orElse(0D), EPSILON);
    }

    private void assertPloidy(double expectedPloidy, boolean alternate, @NotNull final Optional<ModifiableStructuralVariantPloidy> ploidy) {
        assertEquals(expectedPloidy, ploidy.map(ModifiableStructuralVariantPloidy::unweightedImpliedPloidy).orElse(0D), EPSILON);
        if (alternate) {
            assertNotEquals(1d, ploidy.map(ModifiableStructuralVariantPloidy::weight).orElse(0D), EPSILON);
        } else {
            assertEquals(1d, ploidy.map(ModifiableStructuralVariantPloidy::weight).orElse(0D), EPSILON);
        }
    }

    @NotNull
    private PurpleCopyNumber copyNumber(long start, long end, double copyNumber) {
        return PurpleDatamodelTest.createCopyNumber(CHROMOSOME, start, end, copyNumber).build();
    }

    @NotNull
    private ListMultimap<String, PurpleCopyNumber> copyNumbers(PurpleCopyNumber... copyNumbers) {
        final ListMultimap<String, PurpleCopyNumber> result = ArrayListMultimap.create();
        result.putAll(CHROMOSOME, Lists.newArrayList(copyNumbers));
        return result;
    }

}
