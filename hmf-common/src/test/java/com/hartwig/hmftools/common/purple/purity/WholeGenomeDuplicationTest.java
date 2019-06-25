package com.hartwig.hmftools.common.purple.purity;

import static com.hartwig.hmftools.common.purple.purity.WholeGenomeDuplication.MIN_AVERAGE_PLOIDY;
import static com.hartwig.hmftools.common.purple.purity.WholeGenomeDuplication.MIN_DUPLICATED_AUTOSOMES;
import static com.hartwig.hmftools.common.purple.purity.WholeGenomeDuplication.averageMajorAllelePloidy;
import static com.hartwig.hmftools.common.purple.purity.WholeGenomeDuplication.wholeGenomeDuplication;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleDatamodelTest;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WholeGenomeDuplicationTest {

    private static final double EPSILON = 1e-10;

    @Test
    public void testAverageMajorAllelePloidyIsBafCountWeighted() {
        PurpleCopyNumber one = create("1", 1, 1000);
        PurpleCopyNumber two = create("1", 2, 3000);
        PurpleCopyNumber three = create("1", 3, 0);

        assertEquals(1.75, averageMajorAllelePloidy(Lists.newArrayList(one, two, three)), EPSILON);
    }

    @Test
    public void testSexChromosomesExcluded() {
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        for (int i = 1; i < MIN_DUPLICATED_AUTOSOMES; i++) {
            copyNumbers.add(create(String.valueOf(i), MIN_AVERAGE_PLOIDY, 1));
        }
        assertFalse(wholeGenomeDuplication(copyNumbers));

        copyNumbers.add(create(String.valueOf("X"), MIN_AVERAGE_PLOIDY, 1));
        copyNumbers.add(create(String.valueOf("Y"), MIN_AVERAGE_PLOIDY, 1));
        assertFalse(wholeGenomeDuplication(copyNumbers));

        copyNumbers.add(create(String.valueOf(MIN_DUPLICATED_AUTOSOMES), MIN_AVERAGE_PLOIDY, 1));
        assertTrue(wholeGenomeDuplication(copyNumbers));
    }

    @NotNull
    private static PurpleCopyNumber create(String chromosome, double majorAllelePloidy, int bafCount) {
        return PurpleDatamodelTest.createCopyNumber(chromosome, 1, 1, majorAllelePloidy).averageActualBAF(1).bafCount(bafCount).build();
    }
}
