package com.hartwig.hmftools.common.genome.chromosome;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeFunctions;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.purple.Gender;

public enum HumanChromosome implements Chromosome
{
    _1(true),
    _2(true),
    _3(true),
    _4(true),
    _5(true),
    _6(true),
    _7(true),
    _8(true),
    _9(true),
    _10(true),
    _11(true),
    _12(true),
    _13(true),
    _14(true),
    _15(true),
    _16(true),
    _17(true),
    _18(true),
    _19(true),
    _20(true),
    _21(true),
    _22(true),
    _X(false),
    _Y(false);

    private final boolean mIsAutosome;
    private final String mName;

    HumanChromosome(final boolean isAutosome)
    {
        mIsAutosome = isAutosome;
        mName = name().substring(1).intern();
    }

    public static final String CHR_PREFIX = "chr";
    public static final String ENUM_PREFIX = "_";

    @Override
    public boolean isAutosome()
    {
        return mIsAutosome;
    }

    @Override
    public boolean isAllosome() { return !mIsAutosome; }

    public boolean matches(final String chromosome) { return fromString(chromosome) == this; }

    public static Chromosome valueOf(final GenomePosition position)
    {
        return fromString(position.chromosome());
    }

    public static Chromosome valueOf(final GenomeRegion region)
    {
        return fromString(region.chromosome());
    }

    public static HumanChromosome fromString(final String chromosome)
    {
        if(chromosome.toLowerCase().startsWith(CHR_PREFIX))
            return HumanChromosome.valueOf(ENUM_PREFIX + chromosome.substring(3));

        return HumanChromosome.valueOf(ENUM_PREFIX + chromosome);
    }

    public static boolean contains(final String chromosome)
    {
        final String trimmedContig = RefGenomeFunctions.stripChrPrefix(chromosome);
        if(isNumeric(trimmedContig))
        {
            final int integerContig = Integer.parseInt(trimmedContig);
            return integerContig >= 1 && integerContig <= 22;
        }

        return trimmedContig.equals("X") || trimmedContig.equals("Y");
    }

    public int intValue()
    {
        return this.ordinal() + 1;
    }

    public boolean isDiploid(final Gender gender)
    {
        return isAutosome() || (gender != Gender.MALE && this.equals(_X));
    }

    @Override
    public String toString() { return mName; }

    public static boolean isShortArm(final String chromosome)
    {
        return chromosome.equals("13") || chromosome.equals("14") || chromosome.equals("15")
                || chromosome.equals("21") || chromosome.equals("22");
    }

    public boolean isShortArm() { return isShortArm(mName); }

    private static boolean isNumeric(String str)
    {
        for(int i = 0; i < str.length(); i++)
        {
            if(!Character.isDigit(str.charAt(i)))
            {
                return false;
            }
        }
        return true;
    }

    public static boolean lowerChromosome(final String chr, final String otherChr)
    {
        return chromosomeRank(chr) < chromosomeRank(otherChr);
    }

    public static final int INVALID_CHR_RANK = 26;

    public static int chromosomeRank(final String chromosome)
    {
        String chrTrimmed = RefGenomeFunctions.stripChrPrefix(chromosome);

        if(chrTrimmed.equalsIgnoreCase("X"))
        {
            return 23;
        }
        else if(chrTrimmed.equalsIgnoreCase("Y"))
        {
            return 24;
        }
        else if(chrTrimmed.equalsIgnoreCase("MT") || chrTrimmed.equalsIgnoreCase("M"))
        {
            return 25;
        }
        else
        {
            try
            {
                return Integer.parseInt(chrTrimmed);
            }
            catch(NumberFormatException e)
            {
                return INVALID_CHR_RANK;
            }
        }
    }
}
