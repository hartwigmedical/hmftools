package com.hartwig.hmftools.common.purple.gccontent;

import java.util.Map;
import java.util.TreeSet;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

public class ReadMedian {

    private final Map<Integer, ReadCount> readCount;

    public ReadMedian() {
        this.readCount = Maps.newHashMap();
    }

    public void addRead(int read) {
        ReadCount count = readCount.get(read);
        if (count == null) {
            readCount.put(read, new ReadCount(read));
        } else {
            count.increment();
        }
    }

    public int median() {

        TreeSet<ReadCount> sortedSet = Sets.newTreeSet(readCount.values());
        if (sortedSet.isEmpty()) {
            return 0;
        }

        if (sortedSet.size() == 1) {
            return sortedSet.pollFirst().read();
        }

        @NotNull
        ReadCount first = sortedSet.pollFirst();
        int firstCount = first.count();

        @NotNull
        ReadCount last = sortedSet.pollLast();
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

    private static class ReadCount implements Comparable<ReadCount> {
        private final int read;
        private int count;

        private ReadCount(final int read) {
            this.read = read;
            this.count = 1;
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
        public int compareTo(@NotNull final ReadCount other) {
            return Integer.compare(read, other.read);
        }
    }
}
