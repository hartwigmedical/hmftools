package com.hartwig.hmftools.common.purple.ratio;

import java.util.Map;
import java.util.TreeSet;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

class IntegerMedian {

    private final Map<Integer, InnerReadCount> readCount;

    IntegerMedian() {
        this.readCount = Maps.newHashMap();
    }

    void addRead(int read) {
        readCount.computeIfAbsent(read, InnerReadCount::new).increment();
    }

    int median() {

        TreeSet<InnerReadCount> sortedSet = Sets.newTreeSet(readCount.values());
        if (sortedSet.isEmpty()) {
            return 0;
        }

        if (sortedSet.size() == 1) {
            return sortedSet.pollFirst().read();
        }

        @NotNull
        InnerReadCount first = sortedSet.pollFirst();
        int firstCount = first.count();

        @NotNull
        InnerReadCount last = sortedSet.pollLast();
        int lastCount = last.count();

        while (true) {

            if (firstCount > lastCount) {
                if (sortedSet.isEmpty()) {
                    return first.read();
                }
                firstCount -= lastCount;
                last = sortedSet.pollLast();
                lastCount = last.count();
            } else if (lastCount > firstCount) {
                if (sortedSet.isEmpty()) {
                    return last.read();
                }

                lastCount -= firstCount;
                first = sortedSet.pollFirst();
                firstCount = first.count();

            } else {
                if (sortedSet.isEmpty()) {
                    return (first.read() + last.read()) / 2;
                } else if (sortedSet.size() == 1) {
                    return sortedSet.pollFirst().read();
                }

                first = sortedSet.pollFirst();
                firstCount = first.count();
                last = sortedSet.pollLast();
                lastCount = last.count();
            }
        }

    }

    public static class InnerReadCount implements Comparable<InnerReadCount> {
        private final int read;
        private int count;

        private InnerReadCount(final int read) {
            this.read = read;
        }

        public int read() {
            return read;
        }

        public int count() {
            return count;
        }

        void increment() {
            count++;
        }

        @Override
        public int compareTo(@NotNull final InnerReadCount other) {
            return Integer.compare(read, other.read);
        }
    }
}
