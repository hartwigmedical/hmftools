package com.hartwig.hmftools.cobalt.ratio;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.cobalt.targeted.TargetedRatioBuilder;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class TargetedRatioBuilderTest
{
    private static final Chromosome CHROMOSOME = new Chromosome("chr1", 10000);
    private static final double EPSILON = 1e-10;

    @Test
    public void testOnTargetRatio()
    {
        final ListMultimap<Chromosome, ReadRatio> ratios = create(
                readRatio(1001, 0),
                readRatio(2001, 0.5),
                readRatio(11001, 4.0),
                readRatio(12001, 19.5),
                readRatio(23001, 0));

        final Map<GenomePosition, Double> targetEnrichmentRatios = new TreeMap<>();

        targetEnrichmentRatios.put(GenomePositions.create(CHROMOSOME.contig, 2001), 2.0);
        targetEnrichmentRatios.put(GenomePositions.create(CHROMOSOME.contig, 12001), 10.0);
        List<GenomePosition> targetRegions = new ArrayList<>(targetEnrichmentRatios.keySet());

        var ratioBuilder = new TargetedRatioBuilder(targetRegions, targetEnrichmentRatios, ratios);

        ListMultimap<Chromosome, ReadRatio> onTargetRatios = ratioBuilder.onTargetRatios();

        assertEquals(2, onTargetRatios.size());

        ReadRatio readRatio2001 = onTargetRatios.get(CHROMOSOME).get(0);
        ReadRatio readRatio12001 = onTargetRatios.get(CHROMOSOME).get(1);

        assertEquals(2001, readRatio2001.position());

        // ratio = raw ratio / target enrichment / median of raw ratios that overlap with targeted

        // median of the unnormalized gc ratio is 10.0
        // so read ratio = 0.5 / 2.0 / 10 = 0.025
        assertEquals(0.025, readRatio2001.ratio(), EPSILON);

        assertEquals(12001, readRatio12001.position());

        // median of the unnormalized gc ratio is 10.0
        // so read ratio = 19.5 / 10.0 / 10 = 0.195
        assertEquals(0.195, readRatio12001.ratio(), EPSILON);
    }

    @NotNull
    private static ListMultimap<Chromosome, ReadRatio> create(ReadRatio ... readRatios)
    {
        ListMultimap<Chromosome, ReadRatio> ratios = ArrayListMultimap.create();
        ratios.putAll(CHROMOSOME, Arrays.asList(readRatios));
        return ratios;
    }

    @NotNull
    private static ReadRatio readRatio(int position, double ratio)
    {
        return ImmutableReadRatio.builder().chromosome(CHROMOSOME.contig).position(position).ratio(ratio).build();
    }
}
