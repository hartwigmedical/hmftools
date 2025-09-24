package com.hartwig.hmftools.panelbuilder;

import static com.hartwig.hmftools.common.codon.Nucleotides.reverseComplementBases;
import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.DEFAULT_PROBE_QUALITY;
import static com.hartwig.hmftools.panelbuilder.Utils.isDnaSequenceNormal;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.IntStream;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.panelbuilder.probequality.ProbeQualityModel;

import org.jetbrains.annotations.Nullable;

// Handles creation of individual probes.
public class ProbeFactory
{
    private final RefGenomeInterface mRefGenome;
    // Allow these to be null for testing purposes where the quality score doesn't matter.
    @Nullable
    private final ProbeQualityProfile mQualityProfile;
    @Nullable
    private final ProbeQualityModel mQualityModel;

    public ProbeFactory(final RefGenomeInterface refGenome, @Nullable final ProbeQualityProfile qualityProfile,
            @Nullable final ProbeQualityModel qualityModel)
    {
        mRefGenome = refGenome;
        mQualityProfile = qualityProfile;
        mQualityModel = qualityModel;
    }

    // Creates a probe from its region(s)/sequence.
    // Returns empty optional if it's not valid to create a probe at that location.
    public Optional<Probe> createProbe(final SequenceDefinition definition, final TargetMetadata metadata)
    {
        return createProbe(
                definition, metadata,
                () -> buildSequence(definition),
                sequence -> getQualityScore(definition.exactRegionOrNull(), sequence, !definition.isExactRegion()));
    }

    public List<Optional<Probe>> createProbeBatched(final List<SequenceDefinition> definitions, final List<TargetMetadata> metadatas)
    {
        QualityScoreBatcher batcher = new QualityScoreBatcher();
        List<String> sequences = definitions.stream().map(this::buildSequence).toList();
        // TODO: inefficient, don't need to compute quality score for probes which fail validation checks
        sequences.forEach(sequence -> batcher.addQuery(null, sequence, true));
        List<Double> qualityScores = batcher.computeBatch();
        return IntStream.range(0, definitions.size())
                .mapToObj(i -> createProbe(definitions.get(i), metadatas.get(i), () -> sequences.get(i), s -> qualityScores.get(i)))
                .toList();
    }

    private Optional<Probe> createProbe(final SequenceDefinition definition, final TargetMetadata metadata,
            final Supplier<String> getSequence, final Function<String, Double> getQualityScore)
    {
        // Only check properties which are inconvenient for the caller to check in advance.
        // Everything else is expected to be checked by the caller and will generate an exception.

        boolean regionsValid = definition.regions().stream().allMatch(this::isRegionValid);
        if(!regionsValid)
        {
            return Optional.empty();
        }

        String sequence = getSequence.get();
        boolean sequenceValid = sequence.length() == definition.baseLength() && isDnaSequenceNormal(sequence);
        if(!sequenceValid)
        {
            return Optional.empty();
        }

        double qualityScore = getQualityScore.apply(sequence);
        double gcContent = calcGcPercent(sequence);

        return Optional.of(new Probe(definition, sequence, metadata, null, null, qualityScore, gcContent));
    }

    private boolean isRegionValid(final ChrBaseRegion region)
    {
        return region.hasValidPositions() && region.end() <= mRefGenome.getChromosomeLength(region.chromosome());
    }

    private String buildSequence(final SequenceDefinition definition)
    {
        String start = definition.startRegion() == null ? "" : getSequence(definition.startRegion());
        if(definition.startOrientation() == Orientation.REVERSE)
        {
            start = reverseComplementBases(start);
        }
        String insert = definition.insertSequence() == null ? "" : definition.insertSequence();
        String end = definition.endRegion() == null ? "" : getSequence(definition.endRegion());
        if(definition.endOrientation() == Orientation.REVERSE)
        {
            end = reverseComplementBases(end);
        }
        return start + insert + end;
    }

    private String getSequence(final ChrBaseRegion region)
    {
        return mRefGenome.getBaseString(region.chromosome(), region.start(), region.end());
    }

    private double getQualityScore(@Nullable final ChrBaseRegion region, @Nullable final String sequence, boolean allowFallback)
    {
        // TODO: nicer way to do this?
        QualityScoreBatcher batcher = new QualityScoreBatcher();
        batcher.addQuery(region, sequence, allowFallback);
        return batcher.computeBatch().get(0);
    }

    private class QualityScoreBatcher
    {
        private final List<ChrBaseRegion> mRegion = new ArrayList<>();
        private final List<String> mSequence = new ArrayList<>();
        private final List<Boolean> mAllowFallback = new ArrayList<>();

        public void addQuery(@Nullable final ChrBaseRegion region, @Nullable final String sequence, boolean allowFallback)
        {
            mRegion.add(region);
            mSequence.add(sequence);
            mAllowFallback.add(allowFallback);
        }

        public List<Double> computeBatch()
        {
            // Prefer the quality profile first since it's much faster.
            // Fall back to computing the quality score based on the sequence if required.

            ArrayList<OptionalDouble> qualityScores = new ArrayList<>(mRegion.size());
            ArrayList<Integer> indicesForModel = new ArrayList<>();
            ArrayList<byte[]> sequencesForModel = new ArrayList<>();
            for(int i = 0; i < mRegion.size(); ++i)
            {
                ChrBaseRegion region = mRegion.get(i);
                OptionalDouble qualityScore = mQualityProfile != null && region != null
                        ? mQualityProfile.computeQualityScore(region)
                        : OptionalDouble.empty();
                qualityScores.add(qualityScore);

                String sequence = mSequence.get(i);
                if(qualityScore.isEmpty() && mQualityModel != null && sequence != null && mAllowFallback.get(i))
                {
                    indicesForModel.add(i);
                    sequencesForModel.add(sequence.getBytes());
                }
            }

            if(!indicesForModel.isEmpty())
            {
                List<ProbeQualityModel.Result> modelResults = mQualityModel.compute(sequencesForModel);
                for(int i = 0; i < indicesForModel.size(); ++i)
                {
                    qualityScores.set(indicesForModel.get(i), OptionalDouble.of(modelResults.get(i).qualityScore()));
                }
            }

            return qualityScores.stream().map(qualityScore -> qualityScore.orElse(DEFAULT_PROBE_QUALITY)).toList();
        }
    }
}
