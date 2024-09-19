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
    private static final ChromosomalRearrangementsDeterminer CHROMOSOMAL_REARRANGEMENTS_DETERMINER_V37 =
            new ChromosomalRearrangementsDeterminer(OrangeRefGenomeVersion.V37);
    private final RefGenomeCoordinates refCoordinatesV37 = RefGenomeCoordinates.COORDS_37;
    private static final ChromosomalRearrangementsDeterminer CHROMOSOMAL_REARRANGEMENTS_DETERMINER_V38 =
            new ChromosomalRearrangementsDeterminer(OrangeRefGenomeVersion.V38);
    private final RefGenomeCoordinates refCoordinatesV38 = RefGenomeCoordinates.COORDS_38;

    @Test
    public void canDetermine1qTrisomyV37()
    {
        int centromere1 = refCoordinatesV37.centromere("1");
        int endChromosome1 = refCoordinatesV37.length("1");
        canDetermine1qTrisomy(centromere1, endChromosome1, CHROMOSOMAL_REARRANGEMENTS_DETERMINER_V37);
    }

    @Test
    public void canDetermine1qTrisomyV38()
    {
        int centromere1 = refCoordinatesV38.centromere("1");
        int endChromosome1 = refCoordinatesV38.length("1");
        canDetermine1qTrisomy(centromere1, endChromosome1, CHROMOSOMAL_REARRANGEMENTS_DETERMINER_V38);
    }

    @Test
    public void canDetermine1p19qCoDeletionV37()
    {
        int centromere1 = refCoordinatesV37.centromere("1");
        int centromere19 = refCoordinatesV37.centromere("19");
        int endChromosome19 = refCoordinatesV37.length("19");
        canDetermine1p19qCoDeletion(centromere1, centromere19, endChromosome19, CHROMOSOMAL_REARRANGEMENTS_DETERMINER_V37);
    }

    @Test
    public void canDetermine1p19qCoDeletionV38()
    {
        int centromere1 = refCoordinatesV38.centromere("1");
        int centromere19 = refCoordinatesV38.centromere("19");
        int endChromosome19 = refCoordinatesV38.length("19");
        canDetermine1p19qCoDeletion(centromere1, centromere19, endChromosome19, CHROMOSOMAL_REARRANGEMENTS_DETERMINER_V38);
    }

    private void canDetermine1qTrisomy(int centromere1, int endChromosome1,
            ChromosomalRearrangementsDeterminer chromosomalRearrangementsDeterminer)
    {
        int stepSize = (endChromosome1 - centromere1) / 4;
        List<Integer> steps = createSteps(stepSize, centromere1 - 3000000, endChromosome1);

        List<PurpleCopyNumber> arm1Q = Lists.newArrayList(
                copyNumber("1", steps.get(0), steps.get(1), 3),
                copyNumber("1", steps.get(1) + 1, steps.get(2), 2.9),
                copyNumber("1", steps.get(2) + 1, steps.get(3), 2.9)
        );
        List<PurpleCopyNumber> with1qTrisomy = new ArrayList<>(arm1Q);
        with1qTrisomy.add(copyNumber("1", steps.get(3) + 1, endChromosome1, 3.1));
        List<PurpleCopyNumber> without1qTrisomy = new ArrayList<>(arm1Q);
        without1qTrisomy.add(copyNumber("1", steps.get(3) + 1, endChromosome1, 2.5));
        assertTrue(chromosomalRearrangementsDeterminer.determine1qTrisomy(with1qTrisomy));
        assertFalse(chromosomalRearrangementsDeterminer.determine1qTrisomy(without1qTrisomy));
    }

    private void canDetermine1p19qCoDeletion(int centromere1, int centromere19, int endChromosome19,
            ChromosomalRearrangementsDeterminer chromosomalRearrangementsDeterminer)
    {
        int stepSize1 = centromere1 / 4;
        List<Integer> steps1 = createSteps(stepSize1, 1, centromere1 + 3000000);

        int stepSize19 = (endChromosome19 - centromere1) / 4;
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

        List<PurpleCopyNumber> with1pDeletionAndWith19QDeletion = Lists.newArrayList(withDeletion1p);
        with1pDeletionAndWith19QDeletion.addAll(withDeletion19q);
        List<PurpleCopyNumber> with1pDeletionAndWithout19QDeletion = Lists.newArrayList(withDeletion1p);
        with1pDeletionAndWithout19QDeletion.addAll(withoutDeletion19q);
        List<PurpleCopyNumber> without1pDeletionAndWith19QDeletion = Lists.newArrayList(withoutDeletion1p);
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
