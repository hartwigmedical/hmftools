package com.hartwig.hmftools.common.gc;

import java.util.Map;
import java.util.TreeSet;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

class ReadCountMedian {

    private long sum;
    private int count;
    private final Map<Integer, InnerReadCount> readCount;

    ReadCountMedian() {
        this.readCount = Maps.newHashMap();
    }

    void addRead(int read) {
        readCount.computeIfAbsent(read, InnerReadCount::new).increment();
        count++;
        sum += read;
    }

    int mean() {
        return count == 0 ? 0 : (int) (sum / count);
    }

    int median() {
        TreeSet<InnerReadCount> sortedSet = Sets.newTreeSet(readCount.values());
        if (sortedSet.isEmpty()) {
            return 0;
        }

        if (sortedSet.size() == 1) {
            InnerReadCount only = sortedSet.pollFirst();
            assert only != null;
            return only.read();
        }

        InnerReadCount first = sortedSet.pollFirst();
        assert first != null;
        int firstCount = first.count();

        InnerReadCount last = sortedSet.pollLast();
        assert last != null;
        int lastCount = last.count();

        while (true) {
            if (firstCount > lastCount) {
                if (sortedSet.isEmpty()) {
                    return first.read();
                }
                firstCount -= lastCount;
                last = sortedSet.pollLast();
                assert last != null;
                lastCount = last.count();
            } else if (lastCount > firstCount) {
                if (sortedSet.isEmpty()) {
                    return last.read();
                }

                lastCount -= firstCount;
                first = sortedSet.pollFirst();
                assert first != null;
                firstCount = first.count();
            } else {
                if (sortedSet.isEmpty()) {
                    return (first.read() + last.read()) / 2;
                } else if (sortedSet.size() == 1) {
                    InnerReadCount only = sortedSet.pollFirst();
                    assert only != null;
                    return only.read();
                }

                first = sortedSet.pollFirst();
                assert first != null;
                firstCount = first.count();
                last = sortedSet.pollLast();
                assert last != null;
                lastCount = last.count();
            }
        }
    }

    private static class InnerReadCount implements Comparable<InnerReadCount> {
        private final int read;
        private int count;

        private InnerReadCount(final int read) {
            this.read = read;
        }

        int read() {
            return read;
        }

        int count() {
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
