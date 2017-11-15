package com.hartwig.hmftools.common.purple.copynumber.sv;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class StructuralVariantPloidyFactoryTest {

    private static final String CHROMOSOME = "1";

    @Test
    public void testIgnoreInserts() {
        final long svStart = 50;
        final long svEnd = 100;
        final double svVaf = 0.25;

        final ListMultimap<String, PurpleCopyNumber> copyNumbers = create(copyNumber(1, 3333, 4));
        final StructuralVariant insert = sv(svStart, svEnd, StructuralVariantType.INS, svVaf, svVaf);
        assertEquals(0, StructuralVariantPloidyFactory.create(insert, copyNumbers).size());

        final StructuralVariant inv = sv(svStart, svEnd, StructuralVariantType.INV, svVaf, svVaf);
        assertEquals(2, StructuralVariantPloidyFactory.create(inv, copyNumbers).size());

        final StructuralVariant dup = sv(svStart, svEnd, StructuralVariantType.DUP, svVaf, svVaf);
        assertEquals(2, StructuralVariantPloidyFactory.create(dup, copyNumbers).size());

        final StructuralVariant del = sv(svStart, svEnd, StructuralVariantType.DEL, svVaf, svVaf);
        assertEquals(2, StructuralVariantPloidyFactory.create(del, copyNumbers).size());

        final StructuralVariant bnd = sv(svStart, svEnd, StructuralVariantType.BND, svVaf, svVaf);
        assertEquals(2, StructuralVariantPloidyFactory.create(bnd, copyNumbers).size());
    }

    @NotNull
    private StructuralVariant sv(long start, long end, StructuralVariantType type, double startAF, double endAF) {
        return PurpleDatamodelTest.createStructuralVariant(CHROMOSOME, start, CHROMOSOME, end, type).startAF(startAF).endAF(endAF).build();
    }

    @NotNull
    private PurpleCopyNumber copyNumber(long start, long end, double copyNumber) {
        return PurpleDatamodelTest.createCopyNumber(CHROMOSOME, start, end, copyNumber).build();
    }

    @NotNull
    private ListMultimap<String, PurpleCopyNumber> create(PurpleCopyNumber... copyNumbers) {
        final ListMultimap<String, PurpleCopyNumber> result = ArrayListMultimap.create();
        result.putAll(CHROMOSOME, Lists.newArrayList(copyNumbers));
        return result;
    }

    @Ignore
    public void testDel() {
        final StructuralVariant del = PurpleDatamodelTest.createStructuralVariant("1", 3333, "1", 4444, StructuralVariantType.DEL)
                .startAF(0.25)
                .endAF(0.25)
                .build();

        final PurpleCopyNumber first = PurpleDatamodelTest.createCopyNumber("1", 1, 3332, 4).build();
        final PurpleCopyNumber second = PurpleDatamodelTest.createCopyNumber("1", 3333, 4443, 3).build();
        final PurpleCopyNumber third = PurpleDatamodelTest.createCopyNumber("1", 4444, 5000, 4).build();
        final ListMultimap<String, PurpleCopyNumber> copyNumbers = ArrayListMultimap.create();
        copyNumbers.putAll("1", Lists.newArrayList(first, second, third));

        final List<StructuralVariantPloidy> ploidies = StructuralVariantPloidyFactory.create(del, copyNumbers);
        for (StructuralVariantPloidy ploidy : ploidies) {
            System.out.println(ploidy);
        }

    }

    @Ignore
    public void testMultipleDel() {
        final StructuralVariant del1 = PurpleDatamodelTest.createStructuralVariant("1", 3333, "1", 3800, StructuralVariantType.DEL)
                .startAF(0.50)
                .endAF(1d)
                .build();

        final StructuralVariant del2 = PurpleDatamodelTest.createStructuralVariant("1", 3500, "1", 4444, StructuralVariantType.DEL)
                .startAF(1d)
                .endAF(0.5d)
                .build();

        final PurpleCopyNumber first = PurpleDatamodelTest.createCopyNumber("1", 1, 3332, 2).build();
        final PurpleCopyNumber second = PurpleDatamodelTest.createCopyNumber("1", 3333, 3499, 1).build();
        final PurpleCopyNumber third = PurpleDatamodelTest.createCopyNumber("1", 3500, 3799, 0).build();
        final PurpleCopyNumber forth = PurpleDatamodelTest.createCopyNumber("1", 3800, 4443, 1).build();
        final PurpleCopyNumber fifth = PurpleDatamodelTest.createCopyNumber("1", 4444, 5000, 2).build();
        final ListMultimap<String, PurpleCopyNumber> copyNumbers = ArrayListMultimap.create();
        copyNumbers.putAll("1", Lists.newArrayList(first, second, third, forth, fifth));

        final List<StructuralVariantPloidy> ploidies = StructuralVariantPloidyFactory.create(del2, copyNumbers);
        for (StructuralVariantPloidy ploidy : ploidies) {
            System.out.println(ploidy);
        }

    }

}
