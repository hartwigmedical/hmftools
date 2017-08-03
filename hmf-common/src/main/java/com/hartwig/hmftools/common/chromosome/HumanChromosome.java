package com.hartwig.hmftools.common.chromosome;

import java.util.EnumSet;
import java.util.Set;

import com.hartwig.hmftools.common.position.GenomePosition;
import com.hartwig.hmftools.common.purple.gender.Gender;

import org.jetbrains.annotations.NotNull;

public enum HumanChromosome implements Chromosome {
    _1(1, false),
    _2(2, false),
    _3(3, false),
    _4(4, false),
    _5(5, false),
    _6(6, false),
    _7(7, false),
    _8(8, false),
    _9(9, false),
    _10(10, false),
    _11(11, false),
    _12(12, false),
    _13(13, false),
    _14(14, false),
    _15(15, false),
    _16(16, false),
    _17(17, false),
    _18(18, false),
    _19(19, false),
    _20(20, false),
    _21(21, false),
    _22(22, false),
    _X(23, true),
    _Y(24, true),
    _MT(25, false);

    private final int value;
    private final boolean isSexChromosome;

    HumanChromosome(final int value, final boolean isSexChromosome) {
        this.value = value;
        this.isSexChromosome = isSexChromosome;
    }

    @Override
    public boolean isAllosome() {
        return isSexChromosome;
    }

    public static Set<HumanChromosome> autosomes() {
        return EnumSet.complementOf(EnumSet.of(HumanChromosome._X, HumanChromosome._Y, HumanChromosome._MT));
    }

    public static Set<HumanChromosome> allosomes() {
        return EnumSet.of(HumanChromosome._X, HumanChromosome._Y);
    }

    public static Set<HumanChromosome> female() {
        return EnumSet.complementOf(EnumSet.of(HumanChromosome._Y, HumanChromosome._MT));
    }

    public static Chromosome valueOf(@NotNull final GenomePosition position) {
        return fromString(position.chromosome());
    }


    public static Chromosome fromString(@NotNull final String chromosome) {
        return HumanChromosome.valueOf("_" + chromosome);
    }

    public int intValue() {
        return value;
    }

    public boolean isHomologous(@NotNull Gender gender) {
        return gender == Gender.FEMALE || !isAllosome();
    }
}
