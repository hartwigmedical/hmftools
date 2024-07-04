package com.hartwig.hmftools.purple.copynumber;

import static java.util.Collections.singleton;

import static com.hartwig.hmftools.purple.copynumber.StructuralVariantLegsFactoryTest.createLeg;
import static com.hartwig.hmftools.purple.copynumber.StructuralVariantPloidyTest.CHROMOSOME;
import static com.hartwig.hmftools.purple.copynumber.StructuralVariantPloidyTest.PURE;
import static com.hartwig.hmftools.purple.MiscTestUtils.buildPurityAdjuster;

import static org.apache.commons.math3.util.Precision.EPSILON;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.purple.copynumber.sv.ImmutableStructuralVariantLegPloidy;
import com.hartwig.hmftools.purple.copynumber.sv.ImmutableStructuralVariantLegs;
import com.hartwig.hmftools.purple.copynumber.sv.ModifiableStructuralVariantLegPloidy;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidyFactory;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegs;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

@SuppressWarnings("OptionalUsedAsFieldOrParameterType")
public class StructuralVariantPloidyFactoryTest
{
    private static final int AVERAGE_READ_DEPTH = 100;
    private static final double AVERAGE_COPY_NUMBER = 2;

    @NotNull
    public static ImmutableStructuralVariantLegPloidy.Builder svLegPloidy(int orientation, @NotNull final Optional<Double> leftCopyNumber,
            @NotNull final Optional<Double> rightCopyNumber, double ploidy)
    {
        return ImmutableStructuralVariantLegPloidy.builder()
                .chromosome(CHROMOSOME)
                .position(1)
                .orientation((byte) orientation)
                .observedVaf(0.5)
                .adjustedVaf(0.5)
                .alleleFrequency(0.5)
                .homology("")
                .anchoringSupportDistance(0)
                .weight(1)
                .averageImpliedPloidy(ploidy)
                .unweightedImpliedPloidy(ploidy)
                .leftCopyNumber(leftCopyNumber)
                .rightCopyNumber(rightCopyNumber);
    }

    private static final StructuralVariantLegPloidyFactory<PurpleCopyNumber> PURE_PLOIDY_FACTORY = new StructuralVariantLegPloidyFactory<>(
            AVERAGE_READ_DEPTH,
            AVERAGE_COPY_NUMBER,
            PURE,
            PurpleCopyNumber::averageTumorCopyNumber);

    @Test
    public void testSingleValidLeg()
    {
        final StructuralVariantLeg start = createLeg(1000, 1, 0.25);

        final PurpleCopyNumber left = copyNumber(1, 1000, 4);
        final PurpleCopyNumber middle = copyNumber(1001, 2000, 3);
        final PurpleCopyNumber right = copyNumber(2001, 3000, 4);

        final StructuralVariantLegs legs = ImmutableStructuralVariantLegs.builder().start(start).end(Optional.empty()).build();
        final ListMultimap<Chromosome, PurpleCopyNumber> copyNumbers = copyNumbers(left, middle, right);
        final List<StructuralVariantLegPloidy> ploidies = PURE_PLOIDY_FACTORY.create(legs, copyNumbers);
        assertEquals(1, ploidies.size());
        for(StructuralVariantLegPloidy ploidy : ploidies)
        {
            assertEquals(1d, ploidy.averageImpliedPloidy(), EPSILON);
            assertEquals(1d, ploidy.weight(), EPSILON);
        }
    }

    @Test
    public void testBothLegsValid()
    {
        final StructuralVariantLeg start = createLeg(999, 1, 0.25);
        final StructuralVariantLeg end = createLeg(2001, -1, 0.25);

        final PurpleCopyNumber left = copyNumber(1, 1000, 4);
        final PurpleCopyNumber middle = copyNumber(1001, 2000, 3);
        final PurpleCopyNumber right = copyNumber(2001, 3000, 4);

        final StructuralVariantLegs legs = ImmutableStructuralVariantLegs.builder().start(start).end(end).build();
        final ListMultimap<Chromosome, PurpleCopyNumber> copyNumbers = copyNumbers(left, middle, right);
        final List<StructuralVariantLegPloidy> ploidies = PURE_PLOIDY_FACTORY.create(legs, copyNumbers);
        assertEquals(2, ploidies.size());
        for(StructuralVariantLegPloidy ploidy : ploidies)
        {
            assertEquals(1d, ploidy.averageImpliedPloidy(), EPSILON);
            assertEquals(2d, ploidy.weight(), EPSILON);
        }
    }

    @Test
    public void testIncludeNegativeAndZeroCopyNumbersCappedAtZero()
    {
        final StructuralVariantLeg positiveLeg = createLeg(1000, 1, 0.25);
        final StructuralVariantLeg negativeLeg = createLeg(2001, -1, 0.25);
        final PurpleCopyNumber left = copyNumber(1001, 2000, -0.03);
        final PurpleCopyNumber right = copyNumber(2001, 3000, 0);

        ModifiableStructuralVariantLegPloidy leftLegPloidy =
                PURE_PLOIDY_FACTORY.create(positiveLeg, GenomeRegionSelectorFactory.create(singleton(left))).get();
        ModifiableStructuralVariantLegPloidy rightLegPloidy =
                PURE_PLOIDY_FACTORY.create(negativeLeg, GenomeRegionSelectorFactory.create(singleton(right))).get();

        assertFalse(leftLegPloidy.leftCopyNumber().isPresent());
        assertTrue(leftLegPloidy.rightCopyNumber().filter(Doubles::isZero).isPresent());

        assertFalse(rightLegPloidy.leftCopyNumber().isPresent());
        assertTrue(leftLegPloidy.rightCopyNumber().filter(Doubles::isZero).isPresent());
    }

    @Test
    public void testExcludeInfiniteVAF()
    {
        final StructuralVariantLeg leg = createLeg(1001, -1, 1);
        final PurpleCopyNumber left = copyNumber(1, 1000, 3);

        assertFalse(PURE_PLOIDY_FACTORY.create(leg, GenomeRegionSelectorFactory.create(singleton(left))).isPresent());
    }

    @Test
    public void testInferPloidyFromReadDepth()
    {
        int multiplier = 3;

        final StructuralVariantLeg leg = createLeg(1001, -1, 3.5 / 4, AVERAGE_READ_DEPTH * multiplier);
        final PurpleCopyNumber left = copyNumber(1, 1000, 3);
        Optional<ModifiableStructuralVariantLegPloidy> result =
                PURE_PLOIDY_FACTORY.create(leg, GenomeRegionSelectorFactory.create(singleton(left)));
        assertEquals(AVERAGE_COPY_NUMBER * multiplier, result.get().unweightedImpliedPloidy(), EPSILON);
    }

    @Test
    public void testSelectCorrectOrAlternativeCopyNumberForLeg()
    {
        final StructuralVariantLeg positiveLeg = createLeg(1000, 1, 0.25);
        final StructuralVariantLeg negativeLeg = createLeg(2001, -1, 0.25);
        final PurpleCopyNumber left = copyNumber(1, 1000, 4);
        final PurpleCopyNumber middle = copyNumber(1001, 2000, 3);
        final PurpleCopyNumber right = copyNumber(2001, 3000, 4);

        assertPloidy(1, false, PURE_PLOIDY_FACTORY.create(positiveLeg, GenomeRegionSelectorFactory.create(singleton(left))));
        assertPloidy(1, true, PURE_PLOIDY_FACTORY.create(positiveLeg, GenomeRegionSelectorFactory.create(singleton(middle))));

        assertPloidy(1, true, PURE_PLOIDY_FACTORY.create(negativeLeg, GenomeRegionSelectorFactory.create(singleton(middle))));
        assertPloidy(1, false, PURE_PLOIDY_FACTORY.create(negativeLeg, GenomeRegionSelectorFactory.create(singleton(right))));
    }

    @Test
    public void testPurityAdjustedPloidy()
    {
        final StructuralVariantLeg leg = createLeg(1000, 1, 0.5);
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList(copyNumber(1, 1000, 2), copyNumber(1001, 200, 1));

        Optional<ModifiableStructuralVariantLegPloidy> purePloidy =
                PURE_PLOIDY_FACTORY.create(leg, GenomeRegionSelectorFactory.create(copyNumbers));
        assertPloidy(1d, purePloidy);

        final PurityAdjuster diluted = buildPurityAdjuster(Gender.FEMALE, 0.8, 1);
        final StructuralVariantLegPloidyFactory<PurpleCopyNumber> dilutedFactory =
                new StructuralVariantLegPloidyFactory<>(diluted, PurpleCopyNumber::averageTumorCopyNumber);
        Optional<ModifiableStructuralVariantLegPloidy> dilutedPloidy =
                dilutedFactory.create(leg, GenomeRegionSelectorFactory.create(copyNumbers));
        assertPloidy(1.25d, dilutedPloidy);

        final PurityAdjuster male = buildPurityAdjuster(Gender.MALE, 0.8, 1);
        final StructuralVariantLegPloidyFactory<PurpleCopyNumber> maleFactory =
                new StructuralVariantLegPloidyFactory<>(male, PurpleCopyNumber::averageTumorCopyNumber);
        Optional<ModifiableStructuralVariantLegPloidy> malePloidy =
                maleFactory.create(leg, GenomeRegionSelectorFactory.create(copyNumbers));
        assertPloidy(1.125d, malePloidy);
    }

    private void assertPloidy(double expected, @NotNull final Optional<ModifiableStructuralVariantLegPloidy> ploidy)
    {
        assertEquals(expected, ploidy.map(ModifiableStructuralVariantLegPloidy::unweightedImpliedPloidy).orElse(0D), EPSILON);
    }

    private void assertPloidy(double expectedPloidy, boolean alternate,
            @NotNull final Optional<ModifiableStructuralVariantLegPloidy> ploidy)
    {
        assertEquals(expectedPloidy, ploidy.map(ModifiableStructuralVariantLegPloidy::unweightedImpliedPloidy).orElse(0D), EPSILON);
        if(alternate)
        {
            assertNotEquals(1d, ploidy.map(ModifiableStructuralVariantLegPloidy::weight).orElse(0D), EPSILON);
        }
        else
        {
            assertEquals(1d, ploidy.map(ModifiableStructuralVariantLegPloidy::weight).orElse(0D), EPSILON);
        }
    }

    @NotNull
    private static PurpleCopyNumber copyNumber(int start, int end, double copyNumber)
    {
        return PurpleTestUtils.createCopyNumber(CHROMOSOME, start, end, copyNumber).build();
    }

    @NotNull
    private static ListMultimap<Chromosome, PurpleCopyNumber> copyNumbers(@NotNull PurpleCopyNumber... copyNumbers)
    {
        final ListMultimap<Chromosome, PurpleCopyNumber> result = ArrayListMultimap.create();
        result.putAll(HumanChromosome.fromString(CHROMOSOME), Lists.newArrayList(copyNumbers));
        return result;
    }
}
