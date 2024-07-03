package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.purple.PurpleTestUtils.createCopyNumber;
import static com.hartwig.hmftools.purple.fitting.WholeGenomeDuplication.MIN_AVERAGE_PLOIDY;
import static com.hartwig.hmftools.purple.fitting.WholeGenomeDuplication.MIN_DUPLICATED_AUTOSOMES;
import static com.hartwig.hmftools.purple.fitting.WholeGenomeDuplication.averageMajorAlleleCopyNumber;
import static com.hartwig.hmftools.purple.fitting.WholeGenomeDuplication.wholeGenomeDuplication;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.fitting.WholeGenomeDuplication;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class WholeGenomeDuplicationTest
{

    private static final double EPSILON = 1e-10;

    @Test
    public void testKeySet()
    {
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        for(int i = 1; i <= 22; i++)
        {
            copyNumbers.add(create("1", 3, 1000));
            copyNumbers.add(create(String.valueOf(i), 1, 1000));
        }

        assertFalse(WholeGenomeDuplication.wholeGenomeDuplication(copyNumbers));
    }

    @Test
    public void testAverageMajorAlleleCopyNumberIsBafCountWeighted()
    {
        PurpleCopyNumber one = create("1", 1, 1000);
        PurpleCopyNumber two = create("1", 2, 3000);
        PurpleCopyNumber three = create("1", 3, 0);

        assertEquals(1.75, averageMajorAlleleCopyNumber(Lists.newArrayList(one, two, three)), EPSILON);
    }

    @Test
    public void testSexChromosomesExcluded()
    {
        final List<PurpleCopyNumber> copyNumbers = Lists.newArrayList();
        for(int i = 1; i < MIN_DUPLICATED_AUTOSOMES; i++)
        {
            copyNumbers.add(create(String.valueOf(i), MIN_AVERAGE_PLOIDY, 1));
        }
        assertFalse(wholeGenomeDuplication(copyNumbers));

        copyNumbers.add(create("X", MIN_AVERAGE_PLOIDY, 1));
        copyNumbers.add(create("Y", MIN_AVERAGE_PLOIDY, 1));
        assertFalse(wholeGenomeDuplication(copyNumbers));

        copyNumbers.add(create(String.valueOf(MIN_DUPLICATED_AUTOSOMES), MIN_AVERAGE_PLOIDY, 1));
        assertTrue(wholeGenomeDuplication(copyNumbers));
    }

    @NotNull
    private static PurpleCopyNumber create(@NotNull String chromosome, double majorAllelePloidy, int bafCount)
    {
        return createCopyNumber(chromosome, 1, 1, majorAllelePloidy).averageActualBAF(1).bafCount(bafCount).build();
    }
}
