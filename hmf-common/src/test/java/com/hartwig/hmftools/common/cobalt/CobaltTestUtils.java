package com.hartwig.hmftools.common.cobalt;

import static com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes.MIN_Y_COUNT;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;

import org.jetbrains.annotations.NotNull;

public class CobaltTestUtils
{
    @NotNull
    public static CobaltChromosomes female()
    {
        List<MedianRatio> ratios = new ArrayList<>();
        for(int i = 1; i <= 22; i++)
        {
            ratios.add(create(String.valueOf(i), 1, 1));
        }

        ratios.add(create("X", 1, 1));
        return new CobaltChromosomes(ratios);
    }

    @NotNull
    public static CobaltChromosomes male()
    {
        List<MedianRatio> ratios = new ArrayList<>();
        for(int i = 1; i <= 22; i++)
        {
            ratios.add(create(String.valueOf(i), 1, 1));
        }

        ratios.add(create("X", 0.5, 1));
        ratios.add(create("Y", 0.5, MIN_Y_COUNT));
        return new CobaltChromosomes(ratios);
    }

    public static MedianRatio create(final String chromosome, double ratio)
    {
        return create(chromosome, ratio, MIN_Y_COUNT);
    }

    public static MedianRatio create(final String chromosome, double ratio, int count)
    {
        return new MedianRatio(chromosome, ratio, count);
    }

}
