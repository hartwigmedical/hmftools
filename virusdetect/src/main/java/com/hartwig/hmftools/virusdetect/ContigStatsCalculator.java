package com.hartwig.hmftools.virusdetect;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.stream.Collectors.counting;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.partitioningBy;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.virusdetect.VirusConstants.VOTE_CORRECT_BASE_PROBABILITY;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Per-contig support statistics over aligned reads: each read counted once per contig, by its single best alignment
// there. Also attributes each read across the contigs it hits, quantifying strain rivalry within a virus.
public class ContigStatsCalculator
{
    // Held rather than read from the constant directly so tests can inject a controllable value.
    private final double mCorrectBaseProbability;

    public ContigStatsCalculator()
    {
        this(VOTE_CORRECT_BASE_PROBABILITY);
    }

    ContigStatsCalculator(double correctBaseProbability)
    {
        mCorrectBaseProbability = correctBaseProbability;
    }

    public Map<String, ContigStats> compute(String bamFile, ViralReference reference)
    {
        return compute(readAlignments(bamFile), reference);
    }

    public Map<String, ContigStats> compute(Collection<ReadAlignment> alignments, ViralReference reference)
    {
        // Alignments clipping over a contig end are a circular-genome artifact: set them aside (counted per contig) and
        // build the stats from the rest.

        Map<Boolean, List<ReadAlignment>> byOriginClip = alignments.stream()
                .collect(partitioningBy(alignment -> alignment.clipsOverContigEnd(reference.contig(alignment.contig()).length())));
        List<ReadAlignment> filteredAlignments = byOriginClip.get(false);
        List<ReadAlignment> originClipAlignments = byOriginClip.get(true);

        Map<String, ContigAccumulator> accumulators = new HashMap<>();
        filteredAlignments.stream()
                .collect(groupingBy(ReadAlignment::readName))
                .values()
                .forEach(readAlignments -> accumulateRead(readAlignments, accumulators));

        Map<String, Long> originClippedByContig = originClipAlignments.stream().collect(groupingBy(ReadAlignment::contig, counting()));

        return accumulators.entrySet().stream().collect(toMap(
                Map.Entry::getKey, entry -> entry.getValue().toContigStats(
                        entry.getKey(), reference.contig(entry.getKey()).length(),
                        originClippedByContig.getOrDefault(entry.getKey(), 0L).intValue())));
    }

    // Folds one read into the per-contig accumulators: its best alignment (and alignment count) on each contig it hits,
    // plus its cross-contig vote split and strict-win margin.
    private void accumulateRead(List<ReadAlignment> readAlignments, Map<String, ContigAccumulator> accumulators)
    {
        // Collapse BWA -a repeats: per contig keep the best alignment and count how many the read has there.
        Map<String, List<ReadAlignment>> alignmentsByContig = readAlignments.stream().collect(groupingBy(ReadAlignment::contig));

        Map<String, ReadAlignment> bestAlignmentPerContig = alignmentsByContig.entrySet().stream().collect(toMap(
                Map.Entry::getKey, entry -> entry.getValue().stream().reduce(ContigStatsCalculator::chooseBetterAlignment).orElseThrow()));

        alignmentsByContig.forEach((contig, contigAlignments) ->
                {
                    ContigAccumulator accumulator = accumulators.computeIfAbsent(contig, k -> new ContigAccumulator());
                    accumulator.addRead(bestAlignmentPerContig.get(contig), contigAlignments.size());
                }
        );

        addVotes(bestAlignmentPerContig, accumulators);
        addMargin(bestAlignmentPerContig, accumulators);
    }

    // A read's vote splits across the contigs it hits by how well each explains it: every extra base a contig fails to
    // explain multiplies its share by the (pessimistic) chance a base is right, so the closest contig wins most.
    private void addVotes(Map<String, ReadAlignment> bestByContig, Map<String, ContigAccumulator> accumulators)
    {
        // Offsets the exponents for numeric stability.
        int minDivergence = bestByContig.values().stream().mapToInt(ReadAlignment::divergence).min().orElseThrow();

        Map<String, Double> contigWeights = bestByContig.entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> voteWeight(entry.getValue(), minDivergence)));

        double totalWeight = contigWeights.values().stream().mapToDouble(Double::doubleValue).sum();
        contigWeights.forEach((contig, weight) ->
        {
            ContigAccumulator accumulator = accumulators.get(contig);
            accumulator.addVote(weight / totalWeight);
        });
    }

    private double voteWeight(ReadAlignment alignment, int minDivergence)
    {
        return Math.pow(mCorrectBaseProbability, alignment.divergence() - minDivergence);
    }

    // A read aligning to >= 2 contigs is contested. If one strictly beats the rest, its margin (runner-up divergence
    // minus best) is credited to it; a read that ties for best credits no contig, as none stands out.
    private static void addMargin(Map<String, ReadAlignment> bestByContig, Map<String, ContigAccumulator> accumulators)
    {
        if(bestByContig.size() < 2)
        {
            return;
        }

        List<ReadAlignment> ranked = bestByContig.values().stream().sorted(Comparator.comparingInt(ReadAlignment::divergence)).toList();
        ReadAlignment rank0 = ranked.get(0);
        ReadAlignment rank1 = ranked.get(1);

        int margin = rank1.divergence() - rank0.divergence();
        if(margin > 0)
        {
            ContigAccumulator accumulator = accumulators.get(rank0.contig());
            accumulator.addMargin(margin);
        }
    }

    private static List<ReadAlignment> readAlignments(String bamFile)
    {
        List<ReadAlignment> alignments = new ArrayList<>();
        try(SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile)))
        {
            for(SAMRecord record : reader)
            {
                if(!record.getReadUnmappedFlag())
                {
                    alignments.add(ReadAlignment.from(record));
                }
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to read aligned BAM", e);
        }
        return alignments;
    }

    private static ReadAlignment chooseBetterAlignment(ReadAlignment a, ReadAlignment b)
    {
        if(a.alignerScore() != b.alignerScore())
        {
            return a.alignerScore() > b.alignerScore() ? a : b;
        }
        // Deterministic tie-break on location.
        return a.alignmentStart() <= b.alignmentStart() ? a : b;
    }

    // Accumulates one contig's reads (each read's best alignment there), then reduces them to a ContigStats.
    private static class ContigAccumulator
    {
        private final List<ReadAlignment> mAlignments = new ArrayList<>();
        private final List<Integer> mAlignmentCounts = new ArrayList<>();
        private final List<Integer> mMargins = new ArrayList<>();
        private double mVotes;

        private void addRead(ReadAlignment best, int alignmentCount)
        {
            mAlignments.add(best);
            mAlignmentCounts.add(alignmentCount);
        }

        private void addVote(double vote)
        {
            mVotes += vote;
        }

        private void addMargin(int margin)
        {
            mMargins.add(margin);
        }

        private ContigStats toContigStats(String contig, int contigLength, int originClippedReads)
        {
            int[] depth = calculateDepth(contigLength, mAlignments);
            SummaryStats depthSummary = SummaryStats.from(depth);
            int coveredBases = (int) Arrays.stream(depth).filter(d -> d > 0).count();
            int multiAlignReads = (int) mAlignmentCounts.stream().filter(count -> count > 1).count();
            SummaryStats alignerScoreSummary = SummaryStats.from(mAlignments.stream().mapToInt(ReadAlignment::alignerScore).toArray());
            Optional<SummaryStats> marginSummary = mMargins.isEmpty() ? Optional.empty() : Optional.of(SummaryStats.from(mMargins));
            Optional<List<Integer>> marginPercentiles =
                    mMargins.isEmpty() ? Optional.empty() : Optional.of(SummaryStats.percentileCurve(mMargins));
            return new ContigStats(
                    contig, contigLength, mAlignments.size(), multiAlignReads, SummaryStats.from(mAlignmentCounts), originClippedReads,
                    coveredBases, depthSummary, alignerScoreSummary, mVotes, mMargins.size(), marginSummary, marginPercentiles);
        }
    }

    private static int[] calculateDepth(int contigLength, Collection<ReadAlignment> alignments)
    {
        int[] depth = new int[contigLength];
        for(ReadAlignment alignment : alignments)
        {
            for(ReadAlignment.AlignedInterval block : alignment.alignedIntervals())
            {
                int start = max(0, block.referenceStart() - 1);   // intervals are 1-based
                int end = min(contigLength, block.referenceStart() - 1 + block.length());
                for(int position = start; position < end; ++position)
                {
                    ++depth[position];
                }
            }
        }
        return depth;
    }
}
