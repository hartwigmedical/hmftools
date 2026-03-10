package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._10;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._4;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._5;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._6;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._7;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._8;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._9;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;
import static com.hartwig.hmftools.common.purple.GermlineStatus.LIKELY_DIPLOID;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.chromosome.PerChromosomeData;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.purple.region.FittingRegion;

import org.junit.Assert;
import org.junit.Test;

public class AneuploidyDetectorTest
{
    private static final double JUST_ABOVE_CUTOFF = 0.571;
    private static final HumanChromosome[] FIRST_10 = { _1, _2, _3, _4, _5, _6, _7, _8, _9, _10 };

    @Test
    public void regionUsedForAneuploidyDetectionTest()
    {
        PerChromosomeData pcd = new PCD(_1, _2, _3);

        FittingRegion fr1 = fr(_1, 100, 200, GermlineStatus.DIPLOID, 10, 0.5);
        Assert.assertTrue(AneuploidyDetector.useRegionToDetectAneuploidy(pcd, fr1));

        FittingRegion fr2 = fr(_10, 100, 200, GermlineStatus.DIPLOID, 10, 0.5);
        assertFalse(AneuploidyDetector.useRegionToDetectAneuploidy(pcd, fr2));

        FittingRegion fr3 = fr(_1, 100, 200, LIKELY_DIPLOID, 10, 0.5);
        assertFalse(AneuploidyDetector.useRegionToDetectAneuploidy(pcd, fr3));

        FittingRegion fr4 = fr(_1, 100, 200, GermlineStatus.DIPLOID, -1, 0.5);
        assertFalse(AneuploidyDetector.useRegionToDetectAneuploidy(pcd, fr4));

        FittingRegion fr5 = fr(_1, 100, 200, GermlineStatus.DIPLOID, 10, -0.5);
        assertFalse(AneuploidyDetector.useRegionToDetectAneuploidy(pcd, fr5));
    }

    @Test
    public void noRegionsElevatedTest()
    {
        List<FittingRegion> regions = new ArrayList<>();
        for(HumanChromosome chromosome : FIRST_10)
        {
            int position = 1000;
            for(int i = 0; i < 100; i++)
            {
                regions.add(frDiploid(chromosome, position, position + 1000, 10, 0.5));
                position += 2000;
            }
        }
        AneuploidyDetector detector = new AneuploidyDetector(regions, pcdFirst10Chromosomes());
        assertFalse(detector.hasAneuploidy());
    }

    @Test
    public void onePercentOfRegionsElevatedTest()
    {
        List<FittingRegion> regions = new ArrayList<>();
        for(HumanChromosome chromosome : FIRST_10)
        {
            int position = 1000;
            for(int i = 0; i < 100; i++)
            {
                if(i == 10)
                {
                    regions.add(frDiploid(chromosome, position, position + 1000, 10, JUST_ABOVE_CUTOFF));
                }
                else
                {
                    regions.add(frDiploid(chromosome, position, position + 1000, 10, 0.5));
                }
                position += 2000;
            }
        }
        AneuploidyDetector detector = new AneuploidyDetector(regions, pcdFirst10Chromosomes());
        assertTrue(detector.hasAneuploidy());
    }

    @Test
    public void onePercentOfBafCountElevatedTest()
    {
        assertFalse(checkForElevateBafCountOutOfAbout5000(45));
        assertTrue(checkForElevateBafCountOutOfAbout5000(55));
    }

    @Test
    public void elevatedRegionsAreTooShortToCountTest()
    {
        assertFalse(checkFor4PercentElevated(1)); // too short
        assertTrue(checkFor4PercentElevated(2));
    }

    @Test
    public void inputRegionsAreFilteredTest()
    {
        List<FittingRegion> regions = new ArrayList<>();
        for(HumanChromosome chromosome : FIRST_10)
        {
            int position = 1000;
            for(int i = 0; i < 100; i++)
            {
                if(i == 10)
                {
                    regions.add(new FR(chromosome, position, position + 1000, LIKELY_DIPLOID, 10, JUST_ABOVE_CUTOFF, 0.5, 0.5));
                }
                else
                {
                    regions.add(frDiploid(chromosome, position, position + 1000, 10, 0.5));
                }
                position += 2000;
            }
        }
        AneuploidyDetector detector = new AneuploidyDetector(regions, pcdFirst10Chromosomes());
        assertFalse(detector.hasAneuploidy());
    }

    private boolean checkForElevateBafCountOutOfAbout5000(final int elevatedBafPointCount)
    {
        List<FittingRegion> regions = new ArrayList<>();
        for(HumanChromosome chromosome : FIRST_10)
        {
            int position = 1000;
            for(int i = 0; i < 1000; i++)
            {
                if(i == 10)
                {
                    regions.add(frDiploid(chromosome, position, position + 1000, elevatedBafPointCount, JUST_ABOVE_CUTOFF));
                }
                else
                {
                    regions.add(frDiploid(chromosome, position, position + 1000, 5, 0.5));
                }
                position += 2000;
            }
        }
        // total baf count = 999 * 5 + 55 = 5050

        return new AneuploidyDetector(regions, pcdFirst10Chromosomes()).hasAneuploidy();
    }

    private boolean checkFor4PercentElevated(final int multiplier)
    {
        List<FittingRegion> regions = new ArrayList<>();
        for(HumanChromosome chromosome : FIRST_10)
        {
            int position = 1000;
            for(int i = 0; i < 100; i++)
            {
                if(i % 5 == 0)
                {
                    regions.add(frDiploid(chromosome, position, position + 1000, multiplier, JUST_ABOVE_CUTOFF));
                }
                else
                {
                    regions.add(frDiploid(chromosome, position, position + 1000, 5 * multiplier, 0.5));
                }
                position += 2000;
            }
        }
        return new AneuploidyDetector(regions, pcdFirst10Chromosomes()).hasAneuploidy();
    }

    private PerChromosomeData pcdFirst10Chromosomes()
    {
        return new PCD(FIRST_10);
    }

    private FittingRegion fr(HumanChromosome chromosome, int start, int end, GermlineStatus status, int bafCount, double tumorRatio)
    {
        return new FR(chromosome, start, end, status, bafCount, 0.5, 0.5, tumorRatio);
    }

    private FittingRegion frDiploid(HumanChromosome chromosome, int start, int end, int bafCount, double observedBAF)
    {
        return new FR(chromosome, start, end, GermlineStatus.DIPLOID, bafCount, observedBAF, 0.5, 0.5);
    }
}

record FR(HumanChromosome chr, int start, int end,
          GermlineStatus germlineStatus, int bafCount, double observedBAF,
          double observedNormalRatio, double observedTumorRatio) implements FittingRegion
{
    @Override
    public String chromosome()
    {
        return V38.versionedChromosome(chr);
    }
}

class PCD implements PerChromosomeData
{
    private final Set<String> chrNames = new HashSet<>();

    PCD(HumanChromosome... chromosomes)
    {
        for(HumanChromosome chromosome : chromosomes)
        {
            chrNames.add(V38.versionedChromosome(chromosome));
        }
    }

    @Override
    public boolean hasChromosome(final String chromosome)
    {
        return chrNames.contains(chromosome);
    }
}
