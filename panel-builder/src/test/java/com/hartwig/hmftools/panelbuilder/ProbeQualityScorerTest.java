package com.hartwig.hmftools.panelbuilder;

import static java.util.Collections.emptyList;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;
import java.util.OptionalDouble;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.After;
import org.junit.Test;

public class ProbeQualityScorerTest
{
    // Dummy probe data, won't be read by code under test.
    private static final TargetMetadata METADATA = new TargetMetadata(TargetMetadata.Type.CUSTOM, "extra");

    private static final int BATCH_SIZE = 3;

    private final ProbeQualityScorer mScorer = new ProbeQualityScorer(this::computeQualityProfile, this::computeQualityModel, BATCH_SIZE);

    // To test the sequence of calls, populate a list of expected function calls which we pop from and assert on each call.
    private ArrayList<Pair<ChrBaseRegion, OptionalDouble>> mProfileResults = new ArrayList<>();
    private ArrayList<Pair<String, Double>> mModelResults = new ArrayList<>();

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
        assertEquals(emptyList(), mProfileResults);
        assertEquals(emptyList(), mModelResults);
    }

    // Create a probe which will get its quality score from the probe quality profile.
    private static Probe probeForProfile(final ChrBaseRegion region)
    {
        String sequence = "A".repeat(region.baseLength());
        return new Probe(SequenceDefinition.exactRegion(region), sequence, METADATA, null, null, null, 0.0);
    }

    // Create a probe which will get its quality score from the probe quality model.
    private static Probe probeForModel(final String sequence)
    {
        String insertSequence = sequence.substring(1, sequence.length() - 1);
        // Region doesn't really matter as long as the insert sequence is large, to ensure we need to use the probe quality model.
        ChrBaseRegion startRegion = new ChrBaseRegion("1", 1000, 1000);
        ChrBaseRegion endRegion = new ChrBaseRegion("2", 2000, 2000);
        SequenceDefinition definition =
                SequenceDefinition.structuralVariant(startRegion, Orientation.FORWARD, insertSequence, endRegion, Orientation.FORWARD);
        return new Probe(definition, sequence, METADATA, null, null, null, 0.0);
    }

    @Test
    public void testEmptyStream()
    {
        List<Probe> actual = mScorer.computeQualityScores(Stream.empty()).toList();
        List<Probe> expected = emptyList();
        assertEquals(expected, actual);
    }

    @Test
    public void testProfileHasResult()
    {
        // Attempts to get the quality score from the probe quality profile and succeeds.
        ChrBaseRegion region = new ChrBaseRegion("1", 1000, 1120);
        double quality = 0.8;
        Probe probe = probeForProfile(region);
        List<Probe> probes = List.of(probe);
        mProfileResults = new ArrayList<>(List.of(Pair.of(region, OptionalDouble.of(quality))));
        List<Probe> actual = mScorer.computeQualityScores(probes.stream()).toList();
        List<Probe> expected = List.of(probe.withQualityScore(quality));
        assertEquals(expected, actual);
    }

    @Test
    public void testProfileNoResult()
    {
        // Attempts to get the quality score from the probe quality profile and fails, so falls back to the probe quality model.
        ChrBaseRegion region = new ChrBaseRegion("1", 1000, 1120);
        double quality = 0.8;
        Probe probe = probeForProfile(region);
        Stream<Probe> probes = Stream.of(probe);
        mProfileResults = new ArrayList<>(List.of(Pair.of(region, OptionalDouble.empty())));
        mModelResults = new ArrayList<>(List.of(Pair.of(probe.sequence(), quality)));
        List<Probe> actual = mScorer.computeQualityScores(probes).toList();
        List<Probe> expected = List.of(probe.withQualityScore(quality));
        assertEquals(expected, actual);
    }

    @Test
    public void testNeedsModel()
    {
        // Novel sequence requires using the probe quality model.
        String sequence = "ACGTACGT";
        double quality = 0.8;
        Probe probe = probeForModel(sequence);
        Stream<Probe> probes = Stream.of(probe);
        mModelResults = new ArrayList<>(List.of(Pair.of(sequence, quality)));
        List<Probe> actual = mScorer.computeQualityScores(probes).toList();
        List<Probe> expected = List.of(probe.withQualityScore(quality));
        assertEquals(expected, actual);
    }

    @Test
    public void testMultipleFromProfile1()
    {
        // Less than 1 batch of probes.

        ChrBaseRegion region1 = new ChrBaseRegion("1", 1000, 1120);
        double quality1 = 0.8;
        Probe probe1 = probeForProfile(region1);

        ChrBaseRegion region2 = new ChrBaseRegion("2", 2000, 2140);
        double quality2 = 0.9;
        Probe probe2 = probeForProfile(region2);

        Stream<Probe> probes = Stream.of(probe1, probe2);
        mProfileResults = new ArrayList<>(List.of(
                Pair.of(region1, OptionalDouble.of(quality1)),
                Pair.of(region2, OptionalDouble.of(quality2))));
        List<Probe> actual = mScorer.computeQualityScores(probes).toList();
        List<Probe> expected = List.of(probe1.withQualityScore(quality1), probe2.withQualityScore(quality2));
        assertEquals(expected, actual);
    }

    @Test
    public void testMultipleFromProfile2()
    {
        // More than 1 batch of probes.

        ArrayList<Probe> probes = new ArrayList<>();
        ArrayList<Probe> expected = new ArrayList<>();
        for(int i = 0; i < BATCH_SIZE * 10; ++i)
        {
            ChrBaseRegion region = new ChrBaseRegion("1", 1000 * (i + 1), 1000 * (i + 1) + 120 - 1);
            probes.add(probeForProfile(region));
            double quality = 0.1 + i / 1e6;
            mProfileResults.add(Pair.of(region, OptionalDouble.of(quality)));
            expected.add(probes.get(i).withQualityScore(quality));
        }
        List<Probe> actual = mScorer.computeQualityScores(probes.stream()).toList();
        assertEquals(expected, actual);
    }

    @Test
    public void testMultipleFromModel1()
    {
        // Less than 1 batch of probes.

        String sequence1 = "ACGTACGT";
        double quality1 = 0.8;
        Probe probe1 = probeForModel(sequence1);

        String sequence2 = "GTGTGT";
        double quality2 = 0.9;
        Probe probe2 = probeForModel(sequence2);

        Stream<Probe> probes = Stream.of(probe1, probe2);
        mModelResults = new ArrayList<>(List.of(
                Pair.of(sequence1, quality1),
                Pair.of(sequence2, quality2)));
        List<Probe> actual = mScorer.computeQualityScores(probes).toList();
        List<Probe> expected = List.of(probe1.withQualityScore(quality1), probe2.withQualityScore(quality2));
        assertEquals(expected, actual);
    }

    @Test
    public void testMultipleFromModel2()
    {
        // More than 1 batch of probes.

        ArrayList<Probe> probes = new ArrayList<>();
        ArrayList<Probe> expected = new ArrayList<>();
        for(int i = 0; i < BATCH_SIZE * 10; ++i)
        {
            String sequence = "ACG".repeat(i + 1);
            probes.add(probeForModel(sequence));
            double quality = 0.1 + i / 1e6;
            mModelResults.add(Pair.of(sequence, quality));
            expected.add(probes.get(i).withQualityScore(quality));
        }
        List<Probe> actual = mScorer.computeQualityScores(probes.stream()).toList();
        assertEquals(expected, actual);
    }

    // TODO: test mixed data
}
