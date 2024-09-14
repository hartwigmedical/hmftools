package com.hartwig.hmftools.orange.algo.purple;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleTestUtils;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class ChromosomalRearrangementsDeterminerTest
{
    private static final ChromosomalRearrangementsDeterminer chromosomalRearrangementsDeterminer =
            new ChromosomalRearrangementsDeterminer(OrangeRefGenomeVersion.V37);
    private final RefGenomeCoordinates refCoordinates = RefGenomeCoordinates.COORDS_37;

    @Test
    public void canDetermine1qTrisomy()
    {
        int centromere1 = refCoordinates.centromere("1");
        int endChromosome1 = refCoordinates.length("1");

        int stepSize = (endChromosome1 - centromere1) / 4;
        List<Integer> steps = createSteps(stepSize, centromere1 - 3000000, endChromosome1);

        final List<PurpleCopyNumber> arm1Q = Lists.newArrayList(
                copyNumber("1", steps.get(0), steps.get(1), 3),
                copyNumber("1", steps.get(1) + 1, steps.get(2), 2.9),
                copyNumber("1", steps.get(2) + 1, steps.get(3), 2.9)
        );
        final List<PurpleCopyNumber> with1qTrisomy = new ArrayList<>(arm1Q);
        with1qTrisomy.add(copyNumber("1", steps.get(3) + 1, endChromosome1, 3.1));
        final List<PurpleCopyNumber> without1qTrisomy = new ArrayList<>(arm1Q);
        without1qTrisomy.add(copyNumber("1", steps.get(3) + 1, endChromosome1, 2.5));

        assertTrue(chromosomalRearrangementsDeterminer.determine1qTrisomy(with1qTrisomy));
        assertFalse(chromosomalRearrangementsDeterminer.determine1qTrisomy(without1qTrisomy));
    }

    @Test
    public void canDetermine1p19qCoDeletion()
    {
        int centromere1 = refCoordinates.centromere("1");
        int centromere19 = refCoordinates.centromere("19");
        int endChromosome19 = refCoordinates.length("19");

        int stepSize1 = centromere1 / 4;
        List<Integer> steps1 = createSteps(stepSize1, 1, centromere1 + 3000000);

        int stepSize19 = (endChromosome19 - centromere19) / 4;
        List<Integer> steps19 = createSteps(stepSize19, centromere19 - 3000000, endChromosome19);

        List<PurpleCopyNumber> arm1p = Lists.newArrayList(copyNumber("1", 1, steps1.get(1), 0.1), copyNumber("1",
                steps1.get(1) + 1, steps1.get(2), 0.15), copyNumber("1", steps1.get(2) + 1, steps1.get(3), 0.1));
        List<PurpleCopyNumber> arm19q = Lists.newArrayList(copyNumber("19", steps19.get(0), steps19.get(1), 0.1), copyNumber("19",
                steps19.get(1) + 1, steps19.get(2), 0.15), copyNumber("19", steps19.get(2) + 1, steps19.get(3), 0.1));

        List<PurpleCopyNumber> withDeletion1p = new ArrayList<>(arm1p);
        withDeletion1p.add(copyNumber("1", steps1.get(3) + 1, centromere1 + 3000000, 0.1));
        List<PurpleCopyNumber> withoutDeletion1p = new ArrayList<>(arm1p);
        withoutDeletion1p.add(copyNumber("1", steps1.get(3) + 1, centromere1 + 3000000, 0.9));

        List<PurpleCopyNumber> withDeletion19q = new ArrayList<>(arm19q);
        withDeletion19q.add(copyNumber("19", steps19.get(3) + 1, endChromosome19, 0.1));
        List<PurpleCopyNumber> withoutDeletion19q = new ArrayList<>(arm19q);
        withoutDeletion19q.add(copyNumber("19", steps19.get(3) + 1, endChromosome19, 0.9));

        final List<PurpleCopyNumber> with1pDeletionAndWith19QDeletion = Lists.newArrayList(withDeletion1p);
        with1pDeletionAndWith19QDeletion.addAll(withDeletion19q);
        final List<PurpleCopyNumber> with1pDeletionAndWithout19QDeletion = Lists.newArrayList(withDeletion1p);
        with1pDeletionAndWithout19QDeletion.addAll(withoutDeletion19q);
        final List<PurpleCopyNumber> without1pDeletionAndWith19QDeletion = Lists.newArrayList(withoutDeletion1p);
        without1pDeletionAndWith19QDeletion.addAll(withDeletion19q);

        assertTrue(chromosomalRearrangementsDeterminer.determine1p19qCodeletion(with1pDeletionAndWith19QDeletion));
        assertFalse(chromosomalRearrangementsDeterminer.determine1p19qCodeletion(with1pDeletionAndWithout19QDeletion));
        assertFalse(chromosomalRearrangementsDeterminer.determine1p19qCodeletion(without1pDeletionAndWith19QDeletion));
    }

    @NotNull
    private static PurpleCopyNumber copyNumber(String chromosome, int start, int end, double copyNumber)
    {
        return PurpleTestUtils.createCopyNumber(chromosome, start, end, copyNumber).build();
    }

    @NotNull
    private static List<Integer> createSteps(int stepSize, int start, int end)
    {
        List<Integer> steps = new ArrayList<>();
        for(int i = start; i < end; i += stepSize)
        {
            steps.add(i);
        }
        return steps;
    }
}
