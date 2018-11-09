package com.hartwig.hmftools.common.pileup;

import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Pileup implements GenomePosition {

    public abstract int readCount();

    @NotNull
    public abstract String referenceBase();


    @NotNull
    public Set<String> inserts() {
        return countMap().keySet().stream().filter(x -> x.startsWith("+")).map(x -> x.substring(1)).collect(Collectors.toSet());
    }

    @NotNull
    public Set<String> inframeInserts() {
        return inserts().stream().filter(Pileup::inframe).collect(Collectors.toSet());
    }

    @NotNull
    public Set<String> deletes() {
        return countMap().keySet().stream().filter(x -> x.startsWith("-")).map(x -> x.substring(1)).collect(Collectors.toSet());
    }

    @NotNull
    public Set<String> inframeDeletes() {
        return deletes().stream().filter(Pileup::inframe).collect(Collectors.toSet());
    }

    public int mismatchCount(char base) {
        return String.valueOf(base).equals(referenceBase()) ? 0 : countMap().getOrDefault(String.valueOf(base), 0);
    }

    public int mismatchScore(char base) {
        return String.valueOf(base).equals(referenceBase()) ? 0 : scoreMap().getOrDefault(String.valueOf(base), 0);
    }

    public int insertCount(@NotNull final String insert) {
        return countMap().getOrDefault("+" + insert, 0);
    }

    public int insertScore(@NotNull final String insert) {
        return scoreMap().getOrDefault("+" + insert, 0);
    }

    public int deleteCount(@NotNull final String delete) {
        return countMap().getOrDefault("-" + delete, 0);
    }

    public int deleteScore(@NotNull final String delete) {
        return scoreMap().getOrDefault("-" + delete, 0);
    }

    public int referenceCount() {
        return countMap().getOrDefault(String.valueOf(referenceBase()), 0);
    }

    public int referenceScore() {
        return scoreMap().getOrDefault(String.valueOf(referenceBase()), 0);
    }

    public int gMismatchCount() {
        return mismatchCount('G');
    }

    public int aMismatchCount() {
        return mismatchCount('A');
    }

    public int tMismatchCount() {
        return mismatchCount('T');
    }

    public int cMismatchCount() {
        return mismatchCount('C');
    }

    public int insertCount() {
        return countMap().entrySet().stream().filter(x -> x.getKey().startsWith("+")).mapToInt(Map.Entry::getValue).sum();
    }


    public int deleteCount() {
        return countMap().entrySet().stream().filter(x -> x.getKey().startsWith("-")).mapToInt(Map.Entry::getValue).sum();
    }

    public int indelCount() {
        return insertCount() + deleteCount();
    }

    private static boolean inframe(@NotNull final String indel) {
        return (indel.length() - 1) % 3 == 0;
    }

    @NotNull
    protected abstract Map<String, Integer> countMap();

    @NotNull
    protected abstract Map<String, Integer> scoreMap();
}
