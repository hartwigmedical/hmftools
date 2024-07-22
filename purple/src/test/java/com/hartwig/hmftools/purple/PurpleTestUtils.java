package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.REFERENCE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.TUMOR;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.ImmutableCobaltRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.genome.chromosome.CobaltChromosomes;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.segment.Segmentation;

public class PurpleTestUtils
{
    public static final String TUMOR_SAMPLE_ID = "SAMPLE_ID";
    public static final String REF_SAMPLE_ID = "REF_SAMPLE_ID";

    public static PurpleConfig buildDefaultPurpleConfig()
    {
        return buildPurpleConfig(buildDefaultConfigBuilder());
    }

    public static ConfigBuilder buildDefaultConfigBuilder()
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        PurpleConfig.registerConfig(configBuilder);

        configBuilder.setValue(TUMOR, TUMOR_SAMPLE_ID);
        configBuilder.setValue(REFERENCE, REF_SAMPLE_ID);

        return configBuilder;
    }

    public static PurpleConfig buildPurpleConfig(final ConfigBuilder configBuilder)
    {
        return new PurpleConfig("version", configBuilder);
    }

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

    public static CobaltRatio createCobaltRatio(final String chromosome, int position, double tumorRatio)
    {
        return ImmutableCobaltRatio.builder()
                .chromosome(chromosome)
                .position(position)
                .tumorReadDepth(0)
                .referenceReadDepth(0)
                .referenceGCRatio(1)
                .referenceGCDiploidRatio(1)
                .tumorGCRatio(tumorRatio)
                .referenceGcContent(0.5)
                .tumorGcContent(0.5).build();
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


}
