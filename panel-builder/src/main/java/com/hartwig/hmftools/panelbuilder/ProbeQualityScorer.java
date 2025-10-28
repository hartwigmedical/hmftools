package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyIterator;
import static java.util.Objects.requireNonNull;

import static com.hartwig.hmftools.panelbuilder.PanelBuilderConstants.PROBE_QUALITY_PROFILE_MAX_REF_DIFF;
import static com.hartwig.hmftools.panelbuilder.SequenceUtils.sequenceIndelSize;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.OptionalDouble;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.mappability.ProbeQualityProfile;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.panelbuilder.probequality.ProbeQualityModel;

// Helps compute probe quality scores.
// Takes quality scores from the probe quality profile, falling back to the alignment model.
// Also batches computations on the alignment model to make it faster.
public class ProbeQualityScorer
{
    // Use function references rather than exact implementations to allow test mocks.
    private final Function<ChrBaseRegion, OptionalDouble> mComputeQualityProfile;
    private final Function<List<String>, List<Double>> mComputeQualityModel;
    // Aim to batch this many probes together when invoking the probe quality model.
    private final int mModelBatchSize;
    // Buffer at most this many probes total.
    private final int mMaxBufferSize;

    private static final int DEFAULT_MODEL_BATCH_SIZE = 1000;
    private static final int DEFAULT_MAX_BUFFER_SIZE = 10000;

    protected ProbeQualityScorer(final Function<ChrBaseRegion, OptionalDouble> computeQualityProfile,
            final Function<List<String>, List<Double>> computeQualityModel, int modelBatchSize, int maxBufferSize)
    {
        mComputeQualityProfile = computeQualityProfile;
        mComputeQualityModel = computeQualityModel;
        if(modelBatchSize < 1 || maxBufferSize < 1 || modelBatchSize > maxBufferSize)
        {
            throw new IllegalArgumentException("Invalid batching configuration");
        }
        mModelBatchSize = modelBatchSize;
        mMaxBufferSize = maxBufferSize;
    }

    public ProbeQualityScorer(final ProbeQualityProfile qualityProfile, final ProbeQualityModel qualityModel)
    {
        this(
                qualityProfile::computeQualityScore,
                probes -> qualityModel.computeFromSeqString(probes).stream().map(ProbeQualityModel.Result::qualityScore).toList(),
                DEFAULT_MODEL_BATCH_SIZE, DEFAULT_MAX_BUFFER_SIZE);
    }

    public Stream<Probe> computeQualityScores(Stream<Probe> probes)
    {
        return Stream.generate(new QualityScoreGenerator(probes)).takeWhile(Objects::nonNull);
    }

    private class QualityScoreGenerator implements Supplier<Probe>
    {
        private final Iterator<Probe> mSourceIterator;
        private Iterator<Probe> mBatchIterator = emptyIterator();

        public QualityScoreGenerator(Stream<Probe> probes)
        {
            mSourceIterator = probes.iterator();
        }

        @Override
        public Probe get()
        {
            // If we previously processed a batch of probes, then we need to consume that before continuing to maintain probe order.
            if(mBatchIterator.hasNext())
            {
                return mBatchIterator.next();
            }

            // If the next probe's quality score can be computed with the probe quality profile, then yield it immediately.
            // This will be the case most of the time.
            Probe probe;
            if(mSourceIterator.hasNext())
            {
                probe = mSourceIterator.next();
                OptionalDouble qualityScore = tryComputeQualityScoreFromProfile(probe);
                if(qualityScore.isPresent())
                {
                    return probe.withQualityScore(qualityScore.getAsDouble());
                }
            }
            else
            {
                return null;
            }

            // If we got here, then the probe needs the probe quality model, in which case we consume probes until the batch is large enough
            // to process efficiently.
            ArrayList<Probe> buffer = new ArrayList<>(mModelBatchSize);
            buffer.add(probe);
            int modelBatchSize = 1;
            while(mSourceIterator.hasNext() && modelBatchSize < mModelBatchSize && buffer.size() < mMaxBufferSize)
            {
                probe = mSourceIterator.next();
                OptionalDouble qualityScore = tryComputeQualityScoreFromProfile(probe);
                if(qualityScore.isPresent())
                {
                    probe = probe.withQualityScore(qualityScore.getAsDouble());
                }
                else
                {
                    modelBatchSize++;
                }
                buffer.add(probe);
            }
            // Compute the quality scores using the probe quality model and save them for subsequent invocations of the stream.
            List<Probe> resultBatch = computeQualityScoresFromModel(buffer);
            mBatchIterator = resultBatch.iterator();
            return mBatchIterator.next();
        }
    }

    private OptionalDouble tryComputeQualityScoreFromProfile(final Probe probe)
    {
        if(canUseProfile(probe.definition()))
        {
            // Use the worst quality score from the constituent regions. This is most conservative.
            return probe.definition().regions().stream()
                    .map(mComputeQualityProfile)
                    .reduce((q1, q2) ->
                    {
                        if(q1.isPresent() && q2.isPresent())
                        {
                            return OptionalDouble.of(Math.min(q1.getAsDouble(), q2.getAsDouble()));
                        }
                        else
                        {
                            return OptionalDouble.empty();
                        }
                    })
                    .orElseThrow();
        }
        else
        {
            return OptionalDouble.empty();
        }
    }

    private static boolean canUseProfile(final SequenceDefinition sequenceDefinition)
    {
        // If the sequence is very close to the ref genome then there's no need to use the probe quality model.
        // We assume a small perturbation of the ref sequence will not produce a large change in quality score.
        return sequenceIndelSize(sequenceDefinition).orElse(Integer.MAX_VALUE) <= PROBE_QUALITY_PROFILE_MAX_REF_DIFF;
    }

    private List<Probe> computeQualityScoresFromModel(final List<Probe> probes)
    {
        ArrayList<Probe> result = new ArrayList<>(probes.size());
        ArrayList<String> sequences = new ArrayList<>();
        ArrayList<Integer> indices = new ArrayList<>();
        for(int i = 0; i < probes.size(); ++i)
        {
            Probe probe = probes.get(i);
            if(probe.qualityScore() == null)
            {
                sequences.add(requireNonNull(probe.sequence()));
                indices.add(i);
            }
            result.add(probe);
        }

        List<Double> modelResults = mComputeQualityModel.apply(sequences);
        for(int i = 0; i < indices.size(); ++i)
        {
            int index = indices.get(i);
            double qualityScore = modelResults.get(i);
            result.set(index, result.get(index).withQualityScore(qualityScore));
        }

        return result;
    }
}
