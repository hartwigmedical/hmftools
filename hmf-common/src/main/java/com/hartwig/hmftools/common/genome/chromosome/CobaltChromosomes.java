package com.hartwig.hmftools.common.genome.chromosome;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._13;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._15;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._18;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._21;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._X;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._Y;

import java.util.Collection;
import java.util.EnumSet;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public class CobaltChromosomes
{
    private final Gender mGender;
    private final Set<GermlineAberration> mAberrations;
    private final Map<String,CobaltChromosome> mChromosomeMap;

    public static final int MIN_Y_COUNT = 1000;

    public static final double TWO_X_CUTOFF = 0.65;
    public static final double TWO_Y_CUTOFF = 0.75;
    public static final double Y_CUTOFF = 0.05;
    public static final double MOSAIC_X_CUTOFF = 0.8;
    public static final double TRISOMY_CUTOFF = 1.35;
    public static final double TERTRASOMY_CUTOFF = 1.8;

    private static final EnumSet<HumanChromosome> TRISOMY_CHROMOSOMES = EnumSet.of(_13, _15, _18, _21, _X);

    /**
     * Gender:
     * FEMALE := X >= 0.65 && Y < 0.05
     * MALE := !FEMALE
     * Chromosomal Aberrations:
     * MOSAIC_X := FEMALE && X < min(0.8, minAutosomeMedianDepthRatio*)
     * KLINEFELTER (XXY) := MALE && X >= 0.65
     * XYY SYNDROME := MALE && Y >= 0.75
     * TRISOMY_[X,21,13,18,15] := CHR >= 1.4
     * <p>
     * Expected Ratios:
     * FEMALE -> autosome = 1, X = 1, Y = 0
     * MALE -> autosome = 1, allosome = 0.5
     * MOSAIC_X -> X = median X ratio
     * KLINEFELTER (XXY) -> X = 1, Y = 0.5
     * XYY -> X = 0.5, Y = 1
     * TRISOMY_[X,21,13,18,15] -> CHR >= 1.5
     * TERTRASOMY_CUTOFF[9] -> CHR >= 2
     */

    public CobaltChromosomes(final Collection<MedianRatio> unfiltered)
    {
        this(unfiltered, true);
    }

    public CobaltChromosomes(final Collection<MedianRatio> medianRatios, boolean calcAberrations)
    {
        mChromosomeMap = Maps.newHashMap();
        mAberrations = Sets.newHashSet();

        double minAutosomeRatio = 0;
        double yMedian = 0;
        double xMedian = 0;

        for(MedianRatio medianRatio : medianRatios)
        {
            HumanChromosome chromosome = HumanChromosome.fromString(medianRatio.Chromosome);

            if(chromosome.isAutosome())
            {
                minAutosomeRatio = minAutosomeRatio == 0 ? medianRatio.MedianRatio : min(minAutosomeRatio, medianRatio.MedianRatio);
            }
            else if(chromosome == _Y && medianRatio.Count >= MIN_Y_COUNT)
            {
                yMedian = medianRatio.MedianRatio;
            }
            else if(chromosome == _X)
            {
                xMedian = medianRatio.MedianRatio;
            }
        }

        boolean isFemale = Doubles.greaterOrEqual(xMedian, TWO_X_CUTOFF) && Doubles.lessThan(yMedian, Y_CUTOFF);

        mGender = (isFemale) ? Gender.FEMALE : Gender.MALE;

        for(MedianRatio medianRatio : medianRatios)
        {
            HumanChromosome chromosome = HumanChromosome.fromString(medianRatio.Chromosome);

            if(chromosome == _Y && medianRatio.Count < MIN_Y_COUNT)
                continue;

            GermlineAberration aberration = calcAberrations ?
                    aberration(isFemale, chromosome, medianRatio.MedianRatio, minAutosomeRatio) : GermlineAberration.NONE;

            double typicalRatio = typicalRatio(isFemale, chromosome);
            double actualRatio = actualRatio(aberration, typicalRatio, medianRatio.MedianRatio);

            if(aberration != GermlineAberration.NONE)
            {
                mAberrations.add(aberration);
            }

            if(Doubles.positive(typicalRatio))
            {
                CobaltChromosome cobaltChromosome = ImmutableCobaltChromosome.builder()
                        .contig(medianRatio.Chromosome)
                        .typicalRatio(typicalRatio)
                        .actualRatio(actualRatio)
                        .isAllosome(chromosome.isAllosome())
                        .isAutosome(chromosome.isAutosome())
                        .mosiac(aberration == GermlineAberration.MOSAIC_X)
                        .build();

                mChromosomeMap.put(medianRatio.Chromosome, cobaltChromosome);
            }
        }

        if(mAberrations.isEmpty())
        {
            mAberrations.add(GermlineAberration.NONE);
        }
    }

    public boolean hasGermlineAberrations() { return !noAberrations(); }
    public Set<GermlineAberration> germlineAberrations() { return mAberrations; }

    public Gender gender() { return mGender; }

    @NotNull
    public CobaltChromosome get(final String chromosome) { return mChromosomeMap.get(chromosome); }
    public boolean hasChromosome(final String chromosome) { return mChromosomeMap.containsKey(chromosome); }

    public Collection<CobaltChromosome> chromosomes() { return mChromosomeMap.values(); }

    private static double actualRatio(final GermlineAberration aberration, double typicalRatio, double medianRatio)
    {
        switch(aberration)
        {
            case TRISOMY_X:
            case TRISOMY_13:
            case TRISOMY_15:
            case TRISOMY_18:
            case TRISOMY_21:
                return 1.5;
            case TERTRASOMY_9:
                return 2;
            case MOSAIC_X:
                return medianRatio;
            case XYY:
            case KLINEFELTER:
                return 1;
            default:
                return typicalRatio;
        }
    }

    private static GermlineAberration aberration(
            boolean isFemale, final HumanChromosome chromosome, double medianRatio, double minAutosomeRatio)
    {
        if(isTrisomy(chromosome, medianRatio))
        {
            if(chromosome == _13)
                return GermlineAberration.TRISOMY_13;
            else if(chromosome == _15)
                return GermlineAberration.TRISOMY_15;
            else if(chromosome == _18)
                return GermlineAberration.TRISOMY_18;
            else if(chromosome == _21)
                return GermlineAberration.TRISOMY_21;
            else if(chromosome == _X)
                return GermlineAberration.TRISOMY_X;
        }

        if(isXYY(isFemale, chromosome, medianRatio))
            return GermlineAberration.XYY;

        if(isMosiacX(isFemale, chromosome, medianRatio, minAutosomeRatio))
            return GermlineAberration.MOSAIC_X;

        if(isKlinefelterXXY(isFemale, chromosome, medianRatio))
            return GermlineAberration.KLINEFELTER;

        return GermlineAberration.NONE;
    }

    private static boolean isKlinefelterXXY(boolean isFemale, final HumanChromosome chromosome, double medianRatio)
    {
        return !isFemale && chromosome == _X && Doubles.greaterOrEqual(medianRatio, TWO_X_CUTOFF);
    }

    private static boolean isXYY(boolean isFemale, final HumanChromosome chromosome, double medianRatio)
    {
        return !isFemale && chromosome == _Y && Doubles.greaterOrEqual(medianRatio, TWO_Y_CUTOFF);
    }

    private static boolean isMosiacX(boolean isFemale, final HumanChromosome chromosome, double medianRatio, double minAutosomeRatio)
    {
        return isFemale && Doubles.lessThan(medianRatio, min(minAutosomeRatio, MOSAIC_X_CUTOFF)) && chromosome == _X;
    }

    private static boolean isTrisomy(final HumanChromosome chromosome, double medianRatio)
    {
        return Doubles.greaterOrEqual(medianRatio, TRISOMY_CUTOFF) && TRISOMY_CHROMOSOMES.contains(chromosome);
    }

    /*
    private static boolean isTetrasomy(final String chromosome, double medianRatio)
    {
        return isChromosome(chromosome, "9") && Doubles.greaterOrEqual(medianRatio, TERTRASOMY_CUTOFF);
    }
    */

    private static double typicalRatio(boolean isFemale, final HumanChromosome chromosome)
    {
        if(isFemale)
            return chromosome == _Y ? 0d : 1d;

        return chromosome.isAllosome() ? 0.5 : 1d;
    }

    private boolean noAberrations()
    {
        return mAberrations.isEmpty() || (mAberrations.size() == 1 && mAberrations.contains(GermlineAberration.NONE));
    }
}
