package com.hartwig.hmftools.common.chromosome;

import java.util.EnumSet;
import java.util.Set;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public enum HumanChromosome implements Chromosome {
    _1(true, false),
    _2(true, false),
    _3(true, false),
    _4(true, false),
    _5(true, false),
    _6(true, false),
    _7(true, false),
    _8(true, false),
    _9(true, false),
    _10(true, false),
    _11(true, false),
    _12(true, false),
    _13(true, false),
    _14(true, false),
    _15(true, false),
    _16(true, false),
    _17(true, false),
    _18(true, false),
    _19(true, false),
    _20(true, false),
    _21(true, false),
    _22(true, false),
    _X(false, true),
    _Y(false, true),
    _MT(false, false);

    private final boolean isAutosome;
    private final boolean isAllosome;

    HumanChromosome(final boolean isAutosome, boolean isAllosome) {
        this.isAutosome = isAutosome;
        this.isAllosome = isAllosome;
    }

    @Override
    public boolean isAutosome() {
        return isAutosome;
    }

    @Override
    public boolean isAllosome() {
        return isAllosome;
    }

    public static Set<HumanChromosome> female() {
        return EnumSet.complementOf(EnumSet.of(HumanChromosome._Y, HumanChromosome._MT));
    }

    public static Chromosome valueOf(@NotNull final GenomePosition position) {
        return fromString(position.chromosome());
    }

    public static HumanChromosome valueOf(@NotNull final GenomeRegion region) {
        return fromString(region.chromosome());
    }

    public static HumanChromosome fromString(@NotNull final String chromosome) {
        return HumanChromosome.valueOf("_" + chromosome);
    }

    public static boolean contains(@NotNull final String chromosome) {
        try {
            fromString(chromosome);
            return true;
        } catch (IllegalArgumentException e) {
            return false;
        }
    }

    public int intValue() {
        return this.ordinal() + 1;
    }

    public boolean isDiploid(@NotNull Gender gender) {
        return isAutosome() || (gender != Gender.MALE && this.equals(_X));
    }

    @Override
    public String toString() {
        return name().substring(1);
    }
}
