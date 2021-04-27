package com.hartwig.hmftools.common.genome.chromosome;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.Gender;

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
    _Y(false, true);

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

    @NotNull
    public static Chromosome valueOf(@NotNull final GenomePosition position) {
        return fromString(position.chromosome());
    }

    @NotNull
    public static Chromosome valueOf(@NotNull final GenomeRegion region) {
        return fromString(region.chromosome());
    }

    @NotNull
    public static HumanChromosome fromString(@NotNull final String chromosome) {
        if (chromosome.toLowerCase().startsWith("chr")) {
            return HumanChromosome.valueOf("_" + chromosome.substring(3));
        }

        return HumanChromosome.valueOf("_" + chromosome);
    }

    public static boolean contains(@NotNull final String chromosome) {
        final String trimmedConfig = chromosome.toLowerCase().startsWith("chr") ? chromosome.substring(3) : chromosome;
        if (isNumeric(trimmedConfig)) {
            final int integerContig = Integer.parseInt(trimmedConfig);
            return integerContig >= 1 && integerContig <= 22;
        }

        return trimmedConfig.equals("X") || trimmedConfig.equals("Y");
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

    private static boolean isNumeric(String str) {
        for (int i = 0; i < str.length(); i++) {
            if (!Character.isDigit(str.charAt(i))) {
                return false;
            }
        }
        return true;
    }

    public static boolean lowerChromosome(final String chr, final String otherChr)
    {
        return chromosomeRank(chr) < chromosomeRank(otherChr);
    }

    public static int chromosomeRank(final String chromosome)
    {
        if(!HumanChromosome.contains(chromosome))
            return -1;

        if(chromosome.equalsIgnoreCase("X") || chromosome.equalsIgnoreCase("chrX"))
            return 23;
        else if(chromosome.equalsIgnoreCase("Y") || chromosome.equalsIgnoreCase("chrY"))
            return 24;
        else if(chromosome.equalsIgnoreCase("MT") || chromosome.equalsIgnoreCase("chrMT"))
            return 25;
        else
            return Integer.parseInt(chromosome);
    }

}
