package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyIterator;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.OptionalDouble;
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
    private final ProbeQualityProfile mQualityProfile;
    private final ProbeQualityModel mQualityModel;

    // Higher means higher performance efficiency for alignment but more memory usage.
    private static final int BATCH_SIZE = 1000;

    public ProbeQualityScorer(final ProbeQualityProfile qualityProfile, final ProbeQualityModel qualityModel)
    {
        mQualityProfile = qualityProfile;
        mQualityModel = qualityModel;
    }

    public Stream<Probe> getQualityScores(Stream<Probe> probes)
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
            ArrayList<Probe> batch = new ArrayList<>(BATCH_SIZE);
            batch.add(probe);
            while(mSourceIterator.hasNext() && batch.size() < BATCH_SIZE)
            {
                probe = mSourceIterator.next();
                OptionalDouble qualityScore = tryComputeQualityScoreFromProfile(probe);
                if(qualityScore.isPresent())
                {
                    probe = probe.withQualityScore(qualityScore.getAsDouble());
                }
                batch.add(probe);
            }
            // Compute the quality scores using the probe quality model and save them for subsequent invocations of the stream.
            List<Probe> resultBatch = computeQualityScoresFromModel(batch);
            mBatchIterator = resultBatch.iterator();
            return mBatchIterator.next();
        }
    }

    private OptionalDouble tryComputeQualityScoreFromProfile(final Probe probe)
    {
        ChrBaseRegion region = probe.definition().exactRegionOrNull();
        if(region == null)
        {
            return OptionalDouble.empty();
        }
        else
        {
            return mQualityProfile.computeQualityScore(region);
        }
    }

    private List<Probe> computeQualityScoresFromModel(final List<Probe> probes)
    {
        List<Probe> result = new ArrayList<>(probes.size());
        List<byte[]> sequences = new ArrayList<>();
        List<Integer> indices = new ArrayList<>();
        for(int i = 0; i < probes.size(); ++i)
        {
            Probe probe = probes.get(i);
            if(probe.qualityScore() == null)
            {
                sequences.add(probe.sequence().getBytes());
                indices.add(i);
            }
            result.add(probe);
        }

        List<ProbeQualityModel.Result> modelResults = mQualityModel.compute(sequences);
        for(int i = 0; i < modelResults.size(); ++i)
        {
            int index = indices.get(i);
            double qualityScore = modelResults.get(i).qualityScore();
            result.set(index, result.get(index).withQualityScore(qualityScore));
        }

        return result;
    }
}
