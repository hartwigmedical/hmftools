package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyList;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.OptionalDouble;
import java.util.Random;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.MockRefGenome;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.After;
import org.junit.Test;

public class ProbeQualityScorerTest
{
    private static final TargetMetadata METADATA = new TargetMetadata(TargetMetadata.Type.CUSTOM, "extra");

    private static final int BATCH_SIZE = 3;
    private static final int BUFFER_SIZE = 999999;      // I.e. unlimited

    private final ProbeQualityScorer mScorer =
            new ProbeQualityScorer(this::computeQualityProfile, this::computeQualityModel, BATCH_SIZE, BUFFER_SIZE);

    // Helpers for making the test cases more automated.
    private int mProbeCounter = 0;
    private final ArrayList<Probe> mInputProbes = new ArrayList<>();
    private final ArrayList<Probe> mExpectedProbes = new ArrayList<>();

    // To test the sequence of calls, populate a list of expected function calls which we pop from and assert on each call.
    private final ArrayList<Pair<ChrBaseRegion, OptionalDouble>> mProfileResults = new ArrayList<>();
    private final ArrayList<Pair<String, Double>> mModelResults = new ArrayList<>();

    private OptionalDouble computeQualityProfile(final ChrBaseRegion probe)
    {
        Pair<ChrBaseRegion, OptionalDouble> result = mProfileResults.remove(0);
        assertEquals(probe, result.getKey());
        return result.getValue();
    }

    private List<Double> computeQualityModel(final List<String> probes)
    {
        // Assert probes are batched per call as expected.
        assertEquals(Math.min(mModelResults.size(), BATCH_SIZE), probes.size());
        assertTrue(probes.size() < BUFFER_SIZE);

        ArrayList<Double> results = new ArrayList<>();
        for(String probe : probes)
        {
            Pair<String, Double> result = mModelResults.remove(0);
            assertEquals(probe, result.getKey());
            results.add(result.getValue());
        }
        return results;
    }

    @After
    public void after()
    {
        List<Probe> actual = mScorer.computeQualityScores(mInputProbes.stream()).toList();
        assertEquals(mExpectedProbes, actual);
        assertEquals(emptyList(), mProfileResults);
        assertEquals(emptyList(), mModelResults);
    }

    // Probe which attempts to get the quality score from the probe quality profile and succeeds.
    private void probeWithProfileResult()
    {
        mProbeCounter++;
        ChrBaseRegion region = new ChrBaseRegion("1", 1000 * mProbeCounter, 1000 * mProbeCounter + 120 - 1);
        String sequence = "A".repeat(region.baseLength());
        // Probe data other than the region and sequence doesn't really matter.
        Probe probe = new Probe(SequenceDefinition.exactRegion(region), sequence, METADATA, null, null, null, 0.0);
        double quality = 0.1 + mProbeCounter / 1e6;
        mInputProbes.add(probe);
        mProfileResults.add(Pair.of(region, OptionalDouble.of(quality)));
        mExpectedProbes.add(probe.withQualityScore(quality));
    }

    // Probe which attempts to get the quality score from the probe quality profile and fails, so falls back to the probe quality model.
    private void probeWithNoProfileResult()
    {
        mProbeCounter++;
        ChrBaseRegion region = new ChrBaseRegion("1", 1000 * mProbeCounter, 1000 * mProbeCounter + 120 - 1);
        String sequence = "A".repeat(region.baseLength());
        // Probe data other than the region and sequence doesn't really matter.
        Probe probe = new Probe(SequenceDefinition.exactRegion(region), sequence, METADATA, null, null, null, 0.0);
        double quality = 0.1 + mProbeCounter / 1e6;
        mInputProbes.add(probe);
        mProfileResults.add(Pair.of(region, OptionalDouble.empty()));
        mModelResults.add(Pair.of(sequence, quality));
        mExpectedProbes.add(probe.withQualityScore(quality));
    }

    // Probe with a novel sequence that requires using the probe quality model.
    private void probeWithModelResult()
    {
        mProbeCounter++;
        String sequence = MockRefGenome.generateRandomBases(mProbeCounter + 20);
        String insertSequence = sequence.substring(1, sequence.length() - 1);
        // Region doesn't really matter as long as the insert sequence is large, to ensure we need to use the probe quality model.
        ChrBaseRegion startRegion = new ChrBaseRegion("1", 1000 * mProbeCounter, 1000 * mProbeCounter);
        ChrBaseRegion endRegion = new ChrBaseRegion("2", 2000 * mProbeCounter, 2000 * mProbeCounter);
        SequenceDefinition definition =
                SequenceDefinition.structuralVariant(startRegion, Orientation.FORWARD, insertSequence, endRegion, Orientation.FORWARD);
        // Probe data other than the region and sequence doesn't really matter.
        Probe probe = new Probe(definition, sequence, METADATA, null, null, null, 0.0);
        double quality = 0.1 + mProbeCounter / 1e6;
        mInputProbes.add(probe);
        mModelResults.add(Pair.of(sequence, quality));
        mExpectedProbes.add(probe.withQualityScore(quality));
    }

    @Test
    public void testEmptyStream()
    {
        // Nothing to do here, input and output are already initialised to empty.
    }

    @Test
    public void testProfileHasResult()
    {
        // Attempts to get the quality score from the probe quality profile and succeeds.
        probeWithProfileResult();
    }

    @Test
    public void testProfileNoResult()
    {
        // Attempts to get the quality score from the probe quality profile and fails, so falls back to the probe quality model.
        probeWithNoProfileResult();
    }

    @Test
    public void testNeedsModel()
    {
        // Novel sequence requires using the probe quality model.
        probeWithModelResult();
    }

    @Test
    public void testMultipleFromProfile1()
    {
        // Less than 1 batch of probes.
        probeWithProfileResult();
        probeWithProfileResult();
    }

    @Test
    public void testMultipleFromProfile2()
    {
        // More than 1 batch of probes.
        for(int i = 0; i < BATCH_SIZE * 10; ++i)
        {
            probeWithProfileResult();
        }
    }

    @Test
    public void testMultipleFromModel1()
    {
        // Less than 1 batch of probes.
        probeWithModelResult();
        probeWithModelResult();
    }

    @Test
    public void testMultipleFromModel2()
    {
        // More than 1 batch of probes.
        for(int i = 0; i < BATCH_SIZE * 10; ++i)
        {
            probeWithModelResult();
        }
    }

    @Test
    public void testMixed()
    {
        // Interleaved sequence of probes with various quality score calculation pathways.
        probeWithProfileResult();
        probeWithProfileResult();
        probeWithModelResult();
        probeWithModelResult();
        probeWithProfileResult();
        probeWithNoProfileResult();
        probeWithModelResult();
        probeWithNoProfileResult();
        probeWithProfileResult();
        probeWithModelResult();
    }

    @Test
    public void testMixedLarge()
    {
        // Interleaved sequence of probes with various quality score calculation pathways.
        Random random = new Random(42);
        for(int i = 0; i < 10000; ++i)
        {
            int r = random.nextInt(3);
            switch(r)
            {
                case 0:
                    probeWithProfileResult();
                    break;
                case 1:
                    probeWithNoProfileResult();
                    break;
                case 2:
                    probeWithModelResult();
                    break;
                default:
                    throw new RuntimeException();
            }
        }
    }
}
