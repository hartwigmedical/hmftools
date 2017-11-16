package com.hartwig.hmftools.common.purple.copynumber.sv;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.region.GenomeRegionSelector;
import com.hartwig.hmftools.common.region.GenomeRegionSelectorFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class StructuralVariantPloidyFactoryTest {

    private static final String CHROMOSOME = "1";
    private static final double EPSILON = 1e-10;

    @NotNull
    private StructuralVariant sv(long start, long end, StructuralVariantType type, double startAF, double endAF) {
        return PurpleDatamodelTest.createStructuralVariant(CHROMOSOME, start, CHROMOSOME, end, type).startAF(startAF).endAF(endAF).build();
    }

    @NotNull
    private PurpleCopyNumber copyNumber(final String chromosome, long start, long end, double copyNumber) {
        return PurpleDatamodelTest.createCopyNumber(chromosome, start, end, copyNumber).build();
    }

    @NotNull
    private ListMultimap<String, PurpleCopyNumber> create(PurpleCopyNumber... copyNumbers) {
        final ListMultimap<String, PurpleCopyNumber> result = ArrayListMultimap.create();
        result.putAll(CHROMOSOME, Lists.newArrayList(copyNumbers));
        return result;
    }

    private GenomeRegionSelector<PurpleCopyNumber> selector(PurpleCopyNumber... copyNumbers) {
        final List<PurpleCopyNumber> copyNumberList = Lists.newArrayList(copyNumbers);
        return GenomeRegionSelectorFactory.create(copyNumberList);
    }

    @Test
    public void testPurityAdjustedPloidy() {

        final StructuralVariantLeg leg = leg("X", 1001, 1, 0.5);
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList(copyNumber("X", 1, 1000, 2), copyNumber("X", 1001, 200, 1));

        final PurityAdjuster pure = new PurityAdjuster(Gender.FEMALE, 1, 1);
        final StructuralVariantPloidyFactory pureFactory = new StructuralVariantPloidyFactory(pure);
        Optional<ModifiableStructuralVariantPloidy> purePloidy = pureFactory.create(leg, GenomeRegionSelectorFactory.create(copyNumbers));
        assertEquals(1d, purePloidy.map(ModifiableStructuralVariantPloidy::unweightedImpliedPloidy).orElse(0D), EPSILON);

        final PurityAdjuster diluted = new PurityAdjuster(Gender.FEMALE, 0.8, 1);
        final StructuralVariantPloidyFactory dilutedFactory = new StructuralVariantPloidyFactory(diluted);
        Optional<ModifiableStructuralVariantPloidy> dilutedPloidy =
                dilutedFactory.create(leg, GenomeRegionSelectorFactory.create(copyNumbers));
        assertEquals(1.25d, dilutedPloidy.map(ModifiableStructuralVariantPloidy::unweightedImpliedPloidy).orElse(0D), EPSILON);

        final PurityAdjuster male = new PurityAdjuster(Gender.MALE, 0.8, 1);
        final StructuralVariantPloidyFactory maleFactory = new StructuralVariantPloidyFactory(male);
        Optional<ModifiableStructuralVariantPloidy> malePloidy = maleFactory.create(leg, GenomeRegionSelectorFactory.create(copyNumbers));
        assertEquals(1.125d, malePloidy.map(ModifiableStructuralVariantPloidy::unweightedImpliedPloidy).orElse(0D), EPSILON);

    }

    private static StructuralVariantLeg leg(final String chromosome, long position, int orientation, double vaf) {
        return ImmutableStructuralVariantLeg.builder().chromosome(chromosome).position(position).orientation(orientation).vaf(vaf).build();
    }

    //    @Ignore
    //    public void testDel() {
    //        final StructuralVariant del = PurpleDatamodelTest.createStructuralVariant("1", 3333, "1", 4444, StructuralVariantType.DEL)
    //                .startAF(0.25)
    //                .endAF(0.25)
    //                .build();
    //
    //        final PurpleCopyNumber first = PurpleDatamodelTest.createCopyNumber("1", 1, 3332, 4).build();
    //        final PurpleCopyNumber second = PurpleDatamodelTest.createCopyNumber("1", 3333, 4443, 3).build();
    //        final PurpleCopyNumber third = PurpleDatamodelTest.createCopyNumber("1", 4444, 5000, 4).build();
    //        final ListMultimap<String, PurpleCopyNumber> copyNumbers = ArrayListMultimap.create();
    //        copyNumbers.putAll("1", Lists.newArrayList(first, second, third));
    //
    //        final List<StructuralVariantPloidy> ploidies = StructuralVariantPloidyFactory.create(del, copyNumbers);
    //        for (StructuralVariantPloidy ploidy : ploidies) {
    //            System.out.println(ploidy);
    //        }
    //
    //    }

    //    @Ignore
    //    public void testMultipleDel() {
    //        final StructuralVariant del1 = PurpleDatamodelTest.createStructuralVariant("1", 3333, "1", 3800, StructuralVariantType.DEL)
    //                .startAF(0.50)
    //                .endAF(1d)
    //                .build();
    //
    //        final StructuralVariant del2 = PurpleDatamodelTest.createStructuralVariant("1", 3500, "1", 4444, StructuralVariantType.DEL)
    //                .startAF(1d)
    //                .endAF(0.5d)
    //                .build();
    //
    //        final PurpleCopyNumber first = PurpleDatamodelTest.createCopyNumber("1", 1, 3332, 2).build();
    //        final PurpleCopyNumber second = PurpleDatamodelTest.createCopyNumber("1", 3333, 3499, 1).build();
    //        final PurpleCopyNumber third = PurpleDatamodelTest.createCopyNumber("1", 3500, 3799, 0).build();
    //        final PurpleCopyNumber forth = PurpleDatamodelTest.createCopyNumber("1", 3800, 4443, 1).build();
    //        final PurpleCopyNumber fifth = PurpleDatamodelTest.createCopyNumber("1", 4444, 5000, 2).build();
    //        final ListMultimap<String, PurpleCopyNumber> copyNumbers = ArrayListMultimap.create();
    //        copyNumbers.putAll("1", Lists.newArrayList(first, second, third, forth, fifth));
    //
    //        final List<StructuralVariantPloidy> ploidies = StructuralVariantPloidyFactory.create(del2, copyNumbers);
    //        for (StructuralVariantPloidy ploidy : ploidies) {
    //            System.out.println(ploidy);
    //        }
    //
    //    }

}
