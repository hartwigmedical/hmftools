package com.hartwig.hmftools.sage.config;

import static com.hartwig.hmftools.common.utils.ConfigUtils.getConfigValue;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface QualityConfig
{

    String JITTER_PENALTY = "jitter_penalty";
    String JITTER_MIN_REPEAT_COUNT = "jitter_min_repeat_count";
    String BASE_QUAL_FIXED_PENALTY = "base_qual_fixed_penalty";
    String READ_EDGE_FIXED_PENALTY = "read_edge_fixed_penalty";
    String MAP_QUAL_FIXED_PENALTY = "map_qual_fixed_penalty";
    String MAP_QUAL_IMPROPER_PAIR_PENALTY = "map_qual_improper_pair_penalty";
    String MAP_QUAL_READ_EVENTS_PENALTY = "map_qual_read_events_penalty";
    String HIGHLY_POLYMORPHIC_GENES = "highly_polymorphic_genes";

    double DEFAULT_JITTER_PENALTY = 0.25;
    int DEFAULT_JITTER_MIN_REPEAT_COUNT = 3;
    int DEFAULT_BASE_QUAL_FIXED_PENALTY = 12;
    int DEFAULT_READ_EDGE_FIXED_PENALTY = 0;
    int DEFAULT_MAP_QUAL_FIXED_PENALTY = 15;
    int DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY = 15;
    int DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY = 8;
    String DEFAULT_HIGHLY_POLYMORPHIC_GENES = "HLA-A,HLA-B,HLA-C,HLA-DQA1,HLA-DQB1,HLA-DRB1";
    int MAX_HIGHLY_POLYMORPHIC_GENES_QUALITY = 10;

    double jitterPenalty();

    int jitterMinRepeatCount();

    int baseQualityFixedPenalty();

    int distanceFromReadEdgeFixedPenalty();

    int mapQualityFixedPenalty();

    int mapQualityReadEventsPenalty();

    @NotNull
    List<HmfTranscriptRegion> highlyPolymorphicGenes();

    int mapQualityImproperPairPenalty();

    default boolean isHighlyPolymorphic(@NotNull final GenomePosition position)
    {
        for(HmfTranscriptRegion highlyPolymorphicGene : highlyPolymorphicGenes())
        {
            if(highlyPolymorphicGene.contains(position))
            {
                return true;
            }
        }

        return false;
    }

    default int modifiedMapQuality(@NotNull final GenomePosition position, int mapQuality, int readEvents, boolean properPairFlag)
    {
        if(isHighlyPolymorphic(position))
        {
            return Math.min(MAX_HIGHLY_POLYMORPHIC_GENES_QUALITY, mapQuality - mapQualityFixedPenalty());
        }

        int improperPairPenalty = mapQualityImproperPairPenalty() * (properPairFlag ? 0 : 1);
        int distancePenalty = Math.max(0, readEvents - 1) * mapQualityReadEventsPenalty();
        return mapQuality - mapQualityFixedPenalty() - improperPairPenalty - distancePenalty;
    }

    default double modifiedBaseQuality(double baseQuality, int distanceFromReadEdge)
    {
        return Math.min(baseQuality - baseQualityFixedPenalty(), 3 * distanceFromReadEdge - distanceFromReadEdgeFixedPenalty());
    }

    default double jitterPenalty(int repeatCount)
    {
        return (jitterPenalty() * Math.max(0, repeatCount - jitterMinRepeatCount()));
    }

    @NotNull
    static Options createOptions()
    {
        final Options options = new Options();
        options.addOption(HIGHLY_POLYMORPHIC_GENES,
                true,
                "Genes to exclude event distance penalty [" + DEFAULT_HIGHLY_POLYMORPHIC_GENES + "]");

        options.addOption(JITTER_PENALTY,
                true,
                "Penalty to apply to qual score when read context matches with jitter [" + DEFAULT_JITTER_PENALTY + "]");
        options.addOption(JITTER_MIN_REPEAT_COUNT,
                true,
                "Minimum repeat count before applying jitter penalty [" + DEFAULT_JITTER_MIN_REPEAT_COUNT + "]");
        options.addOption(BASE_QUAL_FIXED_PENALTY,
                true,
                "Fixed penalty to apply to base quality [" + DEFAULT_BASE_QUAL_FIXED_PENALTY + "]");
        options.addOption(MAP_QUAL_FIXED_PENALTY, true, "Fixed penalty to apply to map quality [" + DEFAULT_MAP_QUAL_FIXED_PENALTY + "]");
        options.addOption(READ_EDGE_FIXED_PENALTY,
                true,
                "Fixed penalty to apply to distance from read edge [" + DEFAULT_READ_EDGE_FIXED_PENALTY + "]");
        options.addOption(MAP_QUAL_IMPROPER_PAIR_PENALTY,
                true,
                "Penalty to apply to map qual when SAM record does not have the ProperPair flag [" + DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY
                        + "]");
        options.addOption(MAP_QUAL_READ_EVENTS_PENALTY,
                true,
                "Penalty to apply to map qual for additional distance from ref [" + DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY + "]");

        return options;
    }

    @NotNull
    static QualityConfig createConfig(@NotNull final CommandLine cmd, @NotNull final List<HmfTranscriptRegion> allTranscripts)
    {
        final Set<String> highlyPolymorphicGeneNames =
                Arrays.stream(cmd.getOptionValue(HIGHLY_POLYMORPHIC_GENES, DEFAULT_HIGHLY_POLYMORPHIC_GENES).split(","))
                        .collect(Collectors.toSet());

        final List<HmfTranscriptRegion> highlyPolymorphicGeneTranscripts =
                allTranscripts.stream().filter(x -> highlyPolymorphicGeneNames.contains(x.gene())).collect(Collectors.toList());

        return ImmutableQualityConfig.builder()
                .highlyPolymorphicGenes(highlyPolymorphicGeneTranscripts)
                .jitterPenalty(getConfigValue(cmd, JITTER_PENALTY, DEFAULT_JITTER_PENALTY))
                .jitterMinRepeatCount(getConfigValue(cmd, JITTER_MIN_REPEAT_COUNT, DEFAULT_JITTER_MIN_REPEAT_COUNT))
                .baseQualityFixedPenalty(getConfigValue(cmd, BASE_QUAL_FIXED_PENALTY, DEFAULT_BASE_QUAL_FIXED_PENALTY))
                .distanceFromReadEdgeFixedPenalty(getConfigValue(cmd, READ_EDGE_FIXED_PENALTY, DEFAULT_READ_EDGE_FIXED_PENALTY))
                .mapQualityFixedPenalty(getConfigValue(cmd, MAP_QUAL_FIXED_PENALTY, DEFAULT_MAP_QUAL_FIXED_PENALTY))
                .mapQualityReadEventsPenalty(getConfigValue(cmd, MAP_QUAL_READ_EVENTS_PENALTY, DEFAULT_MAP_QUAL_READ_EVENTS_PENALTY))
                .mapQualityImproperPairPenalty(getConfigValue(cmd, MAP_QUAL_IMPROPER_PAIR_PENALTY, DEFAULT_MAP_QUAL_IMPROPER_PAIR_PENALTY))
                .build();
    }
}
