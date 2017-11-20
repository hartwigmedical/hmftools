package com.hartwig.hmftools.common.purple.copynumber.sv;

import static com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantPloidyFactoryTest.copyNumber;
import static com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantPloidyFactoryTest.copyNumbers;
import static com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantPloidyTest.CHROMOSOME;
import static com.hartwig.hmftools.common.purple.copynumber.sv.StructuralVariantPloidyTest.PURE;

import static org.apache.commons.math3.util.Precision.EPSILON;
import static org.junit.Assert.assertEquals;

import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Optional;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantType;

import org.junit.Test;

public class StructuralVariantImpliedCopyNumberTest {

    @Test
    public void testNonsymmetricMultipass() {

        final StructuralVariant firstSV = StructuralVariantLegsFactoryTest.sv(1001, 4001, StructuralVariantType.DEL, 0.25, 0.25);
        final StructuralVariant secondSV = StructuralVariantLegsFactoryTest.sv(2001, 3001, StructuralVariantType.DEL, 1 / 3d, 1 / 3d);

        final PurpleCopyNumber firstCN = copyNumber(1, 1000, 40);
        final PurpleCopyNumber secondCN = copyNumber(1001, 2000, 0);
        final PurpleCopyNumber thirdCN = copyNumber(2001, 3000, 0);
        final PurpleCopyNumber forthCN = copyNumber(3001, 4000, 0);
        final PurpleCopyNumber fifthCN = copyNumber(4001, 5000, 10);

        final List<StructuralVariant> svs = Lists.newArrayList(firstSV, secondSV);
        final ListMultimap<String, PurpleCopyNumber> copyNumbers = copyNumbers(firstCN, secondCN, thirdCN, forthCN, fifthCN);

        final StructuralVariantImpliedCopyNumber victim = new StructuralVariantImpliedCopyNumber(PURE);
        final List<PurpleCopyNumber> result = victim.svImpliedCopyNumber(svs, copyNumbers).get(CHROMOSOME);
        assertEquals(5, result.size());
        assertEquals(40.00, result.get(0).averageTumorCopyNumber(), EPSILON);
        assertEquals(33.75, result.get(1).averageTumorCopyNumber(), EPSILON);
        assertEquals(12.50, result.get(2).averageTumorCopyNumber(), EPSILON);
        assertEquals(03.75, result.get(3).averageTumorCopyNumber(), EPSILON);
        assertEquals(10.00, result.get(4).averageTumorCopyNumber(), EPSILON);
    }

    @Test
    public void testImpliedCopyNumber() {
        final StructuralVariantPloidy left = StructuralVariantPloidyTest.create(1, Optional.of(4d), Optional.empty());
        final StructuralVariantPloidy right = StructuralVariantPloidyTest.create(-1, Optional.empty(), Optional.of(5d));
        final PurpleCopyNumber unknownCopyNumber = PurpleDatamodelTest.createCopyNumber("1", 1, 1000, 0).build();

        final PurpleCopyNumber bothKnown =
                StructuralVariantImpliedCopyNumber.impliedCopyNumber(unknownCopyNumber, Optional.of(left), Optional.of(right));
        assertTrue(bothKnown.inferred());
        assertEquals(3.5, bothKnown.averageTumorCopyNumber(), EPSILON);

        final PurpleCopyNumber leftKnown =
                StructuralVariantImpliedCopyNumber.impliedCopyNumber(unknownCopyNumber, Optional.of(left), Optional.empty());
        assertTrue(bothKnown.inferred());
        assertEquals(3, leftKnown.averageTumorCopyNumber(), EPSILON);

        final PurpleCopyNumber rightKnown =
                StructuralVariantImpliedCopyNumber.impliedCopyNumber(unknownCopyNumber, Optional.empty(), Optional.of(right));
        assertTrue(bothKnown.inferred());
        assertEquals(4, rightKnown.averageTumorCopyNumber(), EPSILON);
    }

}
