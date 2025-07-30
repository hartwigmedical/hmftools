package com.hartwig.hmftools.purple;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.segment.Segmentation;

public class FittingTestUtils
{
    public static CobaltChromosomes buildCobaltChromosomes()
    {
        List<MedianRatio> medianRatios = Lists.newArrayList();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            medianRatios.add(new MedianRatio(chromosome.toString(), 1.0, 1));
        }

        return new CobaltChromosomes(medianRatios);
    }

    public static Segmentation createSegmentation(final ReferenceData referenceData)
    {
        try
        {
            return new Segmentation(referenceData);
        }
        catch(Exception e)
        {
            return null;
        }
    }

    public static final int DEFAULT_REF_DEPTH = 30;
    public static final int DEFAULT_TUMOR_DEPTH = 1000;

    public static AmberBAF createAmberBaf(final String chromosome, final int position, final double tumorBAF, final double normalBAF)
    {
        return new AmberBAF(chromosome, position, tumorBAF, DEFAULT_TUMOR_DEPTH, normalBAF, DEFAULT_REF_DEPTH);
    }

    public static CobaltRatio createCobaltRatio(final Chromosome chromosome, int position, double tumorRatio, double tumorGcContent)
    {
        return new CobaltRatio(
                chromosome.toString(),
                position,
                0,
                0,
                1,
                1,
                0.5,
                tumorRatio,
                tumorGcContent);
    }

    public static ObservedRegion createObservedRegion(final String chromosome, final int start, final int end)
    {
        return new ObservedRegion(
                chromosome, start, end, true, SegmentSupport.NONE, 1, 0.5, 1,
                1, 1, 1, GermlineStatus.DIPLOID, false,
                0.93, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0);
    }

    public static ObservedRegion createObservedRegion(
            final String chromosome, final int start, final int end, double observedBaf, double observedTumorRatio,
            final GermlineStatus germlineStatus, final double tumorCopyNumber)
    {
        return new ObservedRegion(
                chromosome, start, end, true, SegmentSupport.NONE, 10, observedBaf, 10,
                observedTumorRatio, 1, 1, germlineStatus, false,
                0.5, 0, 0, 0, 0, 0,
                0, 2, tumorCopyNumber, observedBaf, tumorCopyNumber, observedBaf);
    }

    public static ObservedRegion createObservedRegion(String chromosome, int start, int end, int bafCount, double observedTumorRatio,
            GermlineStatus germlineStatus)
    {
        return new ObservedRegion(
                chromosome, start, end, true, SegmentSupport.NONE, bafCount, 0.5, 10,
                observedTumorRatio, 1, 1, germlineStatus, false,
                0.5, 0, 0, 0, 0, 0,
                0, 2, 2, 0.5, 0, 0.5);
    }

    public static ObservedRegion createDefaultFittedRegion(final String chromosome, final int start, final int end)
    {
        return new ObservedRegion(
                chromosome, start, end, true, SegmentSupport.NONE, 1, 0.5, 1,
                1, 1, 1, GermlineStatus.DIPLOID, false,
                0.93, 0, 0, 0, 0, 0,
                0, 2, 2, 0.5, 0, 0);
    }

    public static PurityAdjuster buildPurityAdjuster(final Gender gender, final double purity, final double normFactor)
    {
        Map<String, Double> observedRatioMap = Maps.newHashMap();

        for(HumanChromosome chromosome : HumanChromosome.values())
        {
            if(chromosome.isAutosome())
            {
                observedRatioMap.put(chromosome.toString(), 1.0);
            }
            else if(chromosome.equals(HumanChromosome._X))
            {
                if(gender == Gender.MALE)
                {
                    observedRatioMap.put(chromosome.toString(), 0.5);
                }
                else
                {
                    observedRatioMap.put(chromosome.toString(), 1.0);
                }
            }
            else if(chromosome.equals(HumanChromosome._Y))
            {
                if(gender == Gender.MALE)
                {
                    observedRatioMap.put(chromosome.toString(), 0.5);
                }
            }
        }

        return new PurityAdjuster(observedRatioMap, purity, normFactor);
    }
}
