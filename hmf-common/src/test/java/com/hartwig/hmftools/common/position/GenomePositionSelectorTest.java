package com.hartwig.hmftools.common.position;

import static org.junit.Assert.assertEquals;

import java.util.LinkedHashSet;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Optional;
import java.util.function.Consumer;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.chromosome.Chromosome;
import com.hartwig.hmftools.common.collect.Multimaps;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class GenomePositionSelectorTest {

    private GenomePositionSelector<GenomePosition> listSelector;
    private GenomePositionSelector<GenomePosition> iteratorSelector;

    private GenomeRegion region1;
    private GenomeRegion region2;
    private GenomeRegion region3;

    private GenomePosition inRegion1a;
    private GenomePosition inRegion1b;
    private GenomePosition inRegion1c;
    private GenomePosition inRegion2;
    private GenomePosition inRegion3;

    private GenomePosition beforeRegion1;
    private GenomePosition afterRegion2;
    private GenomePosition afterRegion3;

    @Before
    public void setup() {
        region1 = GenomeRegionFactory.create("1", 200, 300);
        region2 = GenomeRegionFactory.create("1", 400, 500);
        region3 = GenomeRegionFactory.create("2", 1, 100);

        inRegion1a = GenomePositions.create("1", 200);
        inRegion1b = GenomePositions.create("1", 250);
        inRegion1c = GenomePositions.create("1", 300);
        inRegion2 = GenomePositions.create("1", 450);
        inRegion3 = GenomePositions.create("2", 50);

        beforeRegion1 = GenomePositions.create("1", 50);
        afterRegion2 = GenomePositions.create("1", 550);
        afterRegion3 = GenomePositions.create("2", 150);

        List<GenomePosition> positions =
                Lists.newArrayList(beforeRegion1, inRegion1a, inRegion1b, inRegion1c, inRegion2, afterRegion2, inRegion3, afterRegion3);

        iteratorSelector = GenomePositionSelectorFactory.create(new LinkedHashSet<>(positions));

        final ListMultimap<Chromosome, GenomePosition> positionMap = Multimaps.fromPositions(positions);
        listSelector = GenomePositionSelectorFactory.create(positionMap);

    }

    @Test
    public void testListSelectorPositionInOrder() {
        assertSelectExisting(listSelector, beforeRegion1);
        assertSelectExisting(listSelector, inRegion1a);
        assertSelectExisting(listSelector, inRegion1b);
        assertSelectExisting(listSelector, inRegion1c);
        assertSelectExisting(listSelector, inRegion2);
        assertSelectExisting(listSelector, afterRegion2);
        assertSelectExisting(listSelector, inRegion3);
        assertSelectExisting(listSelector, afterRegion3);
    }

    @Test
    public void testListSelectorPositionOutOfOrder() {
        assertSelectExisting(listSelector, inRegion1a);
        assertSelectExisting(listSelector, beforeRegion1);
        assertSelectExisting(listSelector, inRegion3);
        assertSelectExisting(listSelector, inRegion1b);
        assertSelectExisting(listSelector, inRegion1c);
        assertSelectExisting(listSelector, inRegion2);
        assertSelectExisting(listSelector, afterRegion2);
        assertSelectExisting(listSelector, afterRegion3);
    }

    @Test
    public void testIteratorSelectorPositionInOrder() {
        assertSelectExisting(iteratorSelector, beforeRegion1);
        assertSelectExisting(iteratorSelector, inRegion1a);
        assertSelectExisting(iteratorSelector, inRegion1b);
        assertSelectExisting(iteratorSelector, inRegion1c);
        assertSelectExisting(iteratorSelector, inRegion2);
        assertSelectExisting(iteratorSelector, afterRegion2);
        assertSelectExisting(iteratorSelector, inRegion3);
        assertSelectExisting(iteratorSelector, afterRegion3);
    }

    @Test(expected = NoSuchElementException.class)
    public void testIteratorSelectorPositionOutOfOrder() {
        assertSelectExisting(iteratorSelector, inRegion1a);
        assertSelectExisting(iteratorSelector, beforeRegion1);
    }

    private void assertSelectExisting(@NotNull final GenomePositionSelector<GenomePosition> victim,
            @NotNull GenomePosition existingPosition) {
        Optional<GenomePosition> result = victim.select(existingPosition);
        assertEquals(result.get(), existingPosition);
    }

    private void assertRegion(@NotNull final GenomePositionSelector<GenomePosition> victim, GenomeRegion region,
            GenomePosition... expectedPositions) {
        final ListConsumer consumer = new ListConsumer();
        victim.select(region, consumer);
        for (int i = 0; i < expectedPositions.length; i++) {
            assertEquals(expectedPositions[i], consumer.positions().get(i));
        }
    }

    @Test
    public void testListSelectInRegion() {
        assertRegion(listSelector, region1, inRegion1a, inRegion1b, inRegion1c);
        assertRegion(listSelector, region2, inRegion2);
        assertRegion(listSelector, region3, inRegion3);
    }

    @Test
    public void testListSelectOutOfOrderRegion() {
        assertRegion(listSelector, region2, inRegion2);
        assertRegion(listSelector, region1, inRegion1a, inRegion1b, inRegion1c);
        assertRegion(listSelector, region1, inRegion1a, inRegion1b, inRegion1c);
        assertRegion(listSelector, region3, inRegion3);
    }

    @Test
    public void testIteratorSelectInRegion() {
        assertRegion(iteratorSelector, region1, inRegion1a, inRegion1b, inRegion1c);
        assertRegion(iteratorSelector, region2, inRegion2);
        assertRegion(iteratorSelector, region3, inRegion3);
    }

    @Test(expected = IllegalArgumentException.class)
    public void testIteratorSelectOutOfOrderRegion() {
        assertRegion(iteratorSelector, region2, inRegion2);
        assertRegion(iteratorSelector, region1, inRegion1a, inRegion1b, inRegion1c);
    }

    @Test(expected = IndexOutOfBoundsException.class)
    public void testIteratorSelectDuplicateRegion() {
        assertRegion(iteratorSelector, region2, inRegion2);
        assertRegion(iteratorSelector, region2, inRegion2);
    }

    private static class ListConsumer implements Consumer<GenomePosition> {

        private final List<GenomePosition> positions = Lists.newArrayList();

        List<GenomePosition> positions() {
            return positions;
        }

        @Override
        public void accept(final GenomePosition genomePosition) {
            positions.add(genomePosition);
        }
    }
}
