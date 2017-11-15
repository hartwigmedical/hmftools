package com.hartwig.hmftools.common.purple.copynumber.sv;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.junit.Test;

public class StructuralVariantCopyNumberTest {

    @Test
    public void testStuff() {
        final StructuralVariant sv1 = PurpleDatamodelTest.createStructuralVariant("11", 38334284, "11", 69626993, StructuralVariantType.INV)
                .startAF(0.8610)
                .endAF(0.884)
                .build();
        final StructuralVariant sv2 = PurpleDatamodelTest.createStructuralVariant("11", 68080885, "11", 69626994, StructuralVariantType.DUP)
                .startAF(0.812)
                .endAF(0.823)
                .build();

        final List<StructuralVariant> structuralVariants = Lists.newArrayList(sv1, sv2);
        //        final PurpleCopyNumber negone = PurpleDatamodelTest.createCopyNumber("11", 38333979, 38334283, 196.136).build();
        final PurpleCopyNumber zero = PurpleDatamodelTest.createCopyNumber("11", 38334284, 38347000, 57.97).build();
        final PurpleCopyNumber first = PurpleDatamodelTest.createCopyNumber("11", 69351788, 69626992, 50.62).build();
        final PurpleCopyNumber second = PurpleDatamodelTest.createCopyNumber("11", 69626993, 69626993, 0).build();
        final PurpleCopyNumber third = PurpleDatamodelTest.createCopyNumber("11", 69626994, 76011780, 2.52).build();
        final ListMultimap<String, PurpleCopyNumber> copyNumbers = ArrayListMultimap.create();
        copyNumbers.putAll("11", Lists.newArrayList(zero, first, second, third));

        final ListMultimap<String, PurpleCopyNumber> result = new StructuralVariantCopyNumber(structuralVariants).calculateSVCopyNumber(copyNumbers);

        List<StructuralVariantPloidy> ploidies = StructuralVariantPloidyFactory.create(structuralVariants, copyNumbers);
        for (StructuralVariantPloidy ploidy : ploidies) {
            System.out.println(ploidy);
            //        }

            for (PurpleCopyNumber purpleCopyNumber : result.values()) {
                System.out.println(purpleCopyNumber);
            }
        }
    }

}
