package com.hartwig.hmftools.common.segmentation;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class DoubleArrayTest extends SegmentationTestBase
{

    @Test
    public void powerSetTest() {
        Set<Integer> set1 = new HashSet<>(List.of(1));
        Set<Set<Integer>> powerSet1 = powerSet(set1);
        assertEquals(2, powerSet1.size());
        assertTrue(powerSet1.contains(new HashSet<>()));
        assertTrue(powerSet1.contains(new HashSet<>(Arrays.asList(1))));

        Set<Integer> set2 = new HashSet<>(Arrays.asList(1, 2));
        Set<Set<Integer>> powerSet2 = powerSet(set2);
        assertEquals(4, powerSet2.size());
        assertTrue(powerSet2.contains(new HashSet<>()));
        assertTrue(powerSet2.contains(new HashSet<>(Arrays.asList(1))));
        assertTrue(powerSet2.contains(new HashSet<>(Arrays.asList(2))));
        assertTrue(powerSet2.contains(new HashSet<>(Arrays.asList(1, 2))));

        Set<Integer> set3 = new HashSet<>(Arrays.asList(1, 2, 3));
        Set<Set<Integer>> powerSet3 = powerSet(set3);
        assertEquals(8, powerSet3.size());
        assertTrue(powerSet3.contains(new HashSet<>()));
        assertTrue(powerSet3.contains(new HashSet<>(Arrays.asList(1))));
        assertTrue(powerSet3.contains(new HashSet<>(Arrays.asList(2))));
        assertTrue(powerSet3.contains(new HashSet<>(Arrays.asList(3))));
        assertTrue(powerSet3.contains(new HashSet<>(Arrays.asList(1, 2))));
        assertTrue(powerSet3.contains(new HashSet<>(Arrays.asList(1, 3))));
        assertTrue(powerSet3.contains(new HashSet<>(Arrays.asList(2, 3))));
        assertTrue(powerSet3.contains(new HashSet<>(Arrays.asList(1, 2, 3))));
    }

    @Test
    public void splitByIndexesTest() {
        List<double[]> result1 = splitByIndexes(d(1, 2), new HashSet<>());
        assertEquals(1, result1.size());
        assertArrayEquals(d(1, 2), result1.get(0), 1e-10);

        List<double[]> result2 = splitByIndexes(d(1, 2), new HashSet<>(Arrays.asList(1)));
        assertEquals(2, result2.size());
        assertArrayEquals(d(1), result2.get(0), 1e-10);
        assertArrayEquals(d(2), result2.get(1), 1e-10);

        List<double[]> result3 = splitByIndexes(d(1, 2, 3), new HashSet<>(Arrays.asList(1)));
        assertEquals(2, result3.size());
        assertArrayEquals(d(1), result3.get(0), 1e-10);
        assertArrayEquals(d(2, 3), result3.get(1), 1e-10);

        List<double[]> result4 = splitByIndexes(d(1, 2, 1), new HashSet<>(Arrays.asList(1)));
        assertEquals(2, result4.size());
        assertArrayEquals(d(1), result4.get(0), 1e-10);
        assertArrayEquals(d(2, 1), result4.get(1), 1e-10);

        List<double[]> result5 = splitByIndexes(d(1, 2, 1), new HashSet<>(Arrays.asList(1, 2)));
        assertEquals(3, result5.size());
        assertArrayEquals(d(1), result5.get(0), 1e-10);
        assertArrayEquals(d(2), result5.get(1), 1e-10);
        assertArrayEquals(d(1), result5.get(2), 1e-10);

        List<double[]> result6 = splitByIndexes(d(1, 2, 3, 5, 7), new HashSet<>(Arrays.asList(2)));
        assertEquals(2, result6.size());
        assertArrayEquals(d(1, 2), result6.get(0), 1e-10);
        assertArrayEquals(d(3, 5, 7), result6.get(1), 1e-10);

        List<double[]> result7 = splitByIndexes(d(1, 2, 3, 5, 7), new HashSet<>(Arrays.asList(2, 4)));
        assertEquals(3, result7.size());
        assertArrayEquals(d(1, 2), result7.get(0), 1e-10);
        assertArrayEquals(d(3, 5), result7.get(1), 1e-10);
        assertArrayEquals(d(7), result7.get(2), 1e-10);
    }

    @Test
    public void splitAtTest() {
        Pair<double[], double[]> result1 = splitAt(d(1, 2), 1);
        assertArrayEquals(d(1), result1.getFirst(), 1e-10);
        assertArrayEquals(d(2), result1.getSecond(), 1e-10);

        Pair<double[], double[]> result2 = splitAt(d(1, 2, 3), 1);
        assertArrayEquals(d(1), result2.getFirst(), 1e-10);
        assertArrayEquals(d(2, 3), result2.getSecond(), 1e-10);

        Pair<double[], double[]> result3 = splitAt(d(1, 2, 3), 2);
        assertArrayEquals(d(1, 2), result3.getFirst(), 1e-10);
        assertArrayEquals(d(3), result3.getSecond(), 1e-10);

        Pair<double[], double[]> result4 = splitAt(d(1, 2, 3, 4, 5), 3);
        assertArrayEquals(d(1, 2, 3), result4.getFirst(), 1e-10);
        assertArrayEquals(d(4, 5), result4.getSecond(), 1e-10);
    }

    @Test
    public void medianAbsoluteDeviationTest() {
        assertEquals(0, medianAbsoluteDeviation(d(1, 2, 1)), 1e-10);
        assertEquals(1.4826, medianAbsoluteDeviation(d(10, 12, 11, 14)), 1e-10);
    }

    /**
     * Computes the power set of a set.
     */
    private <T> Set<Set<T>> powerSet(Set<T> set) {
        Set<Set<T>> result = new HashSet<>();
        result.add(new HashSet<>());
        for (T element : set) {
            Set<T> withoutElement = new HashSet<>(set);
            withoutElement.remove(element);
            Set<Set<T>> power = powerSet(withoutElement);
            for (Set<T> subset : power) {
                result.add(subset);
                Set<T> withElement = new HashSet<>(subset);
                withElement.add(element);
                result.add(withElement);
            }
        }
        return result;
    }

    private List<double[]> splitByIndexes(double[] array, Set<Integer> indexes) {
        for (int index : indexes) {
            if (index <= 0 || index >= array.length) {
                throw new IllegalArgumentException("Index must be > 0 and < array.length");
            }
        }
        if (indexes.isEmpty()) {
            return Arrays.asList(array);
        }

        int greatestIndex = Integer.MIN_VALUE;
        for (int index : indexes) {
            if (index > greatestIndex) {
                greatestIndex = index;
            }
        }
        Pair<double[], double[]> headTail = splitAt(array, greatestIndex);
        Set<Integer> remainingIndexes = new HashSet<>(indexes);
        remainingIndexes.remove(greatestIndex);
        List<double[]> headSplits = splitByIndexes(headTail.getFirst(), remainingIndexes);
        List<double[]> result = new ArrayList<>(headSplits);
        result.add(headTail.getSecond());
        return result;
    }

    private Pair<double[], double[]> splitAt(double[] array, int index) {
        if (index <= 0 || index >= array.length) {
            throw new IllegalArgumentException("Index must be > 0 and < array.length");
        }
        double[] left = new double[index];
        double[] right = new double[array.length - index];
        for (int i = 0; i < array.length; i++) {
            if (i < index) {
                left[i] = array[i];
            } else {
                right[i - index] = array[i];
            }
        }
        return new Pair<>(left, right);
    }

    private double medianAbsoluteDeviation(double[] array) {
        return Stats.medianAbsoluteDeviation(array);
    }

    private static class Pair<F, S> {
        private final F first;
        private final S second;

        public Pair(F first, S second) {
            this.first = first;
            this.second = second;
        }

        public F getFirst() {
            return first;
        }

        public S getSecond() {
            return second;
        }
    }
}
