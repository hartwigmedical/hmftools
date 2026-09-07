package com.hartwig.hmftools.virusdetect;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.virusdetect.VirusConstants.VOTE_CORRECT_BASE_PROBABILITY;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Per-contig support statistics over an aligned BAM: each read counted once per contig, by its single best alignment
// there. Also attributes each read across the contigs it hits, quantifying strain rivalry within a virus.
public class ContigStatsCalculator
{
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
        Map<String, Map<String, SAMRecord>> alignmentsByRead = readBestAlignments(bamFile);

        Map<String, List<SAMRecord>> recordsByContig = new HashMap<>();
        Map<String, Double> votesByContig = new HashMap<>();
        Map<String, List<Integer>> marginsByContig = new HashMap<>();

        for(Map<String, SAMRecord> byContig : alignmentsByRead.values())
        {
            byContig.forEach((contig, record) -> recordsByContig.computeIfAbsent(contig, k -> new ArrayList<>()).add(record));
            addVotes(byContig, votesByContig);
            addMargin(byContig, marginsByContig);
        }

        Map<String, ContigStats> stats = new HashMap<>();
        for(Map.Entry<String, List<SAMRecord>> entry : recordsByContig.entrySet())
        {
            String contig = entry.getKey();
            stats.put(contig, computeContig(
                    contig, reference.contig(contig).length(), entry.getValue(),
                    votesByContig.getOrDefault(contig, 0.0), marginsByContig.get(contig)));
        }
        return stats;
    }

    // Best (highest-scoring) alignment of each read on each contig; a read holds at most one record per contig.
    private static Map<String, Map<String, SAMRecord>> readBestAlignments(String bamFile)
    {
        Map<String, Map<String, SAMRecord>> alignmentsByRead = new HashMap<>();
        try(SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile)))
        {
            for(SAMRecord record : reader)
            {
                if(record.getReadUnmappedFlag())
                {
                    continue;
                }
                alignmentsByRead.computeIfAbsent(record.getReadName(), k -> new HashMap<>())
                        .merge(record.getReferenceName(), record, ContigStatsCalculator::better);
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to read aligned BAM", e);
        }
        return alignmentsByRead;
    }

    // A read's vote splits across the contigs it hits by how well each explains it: every extra base a contig fails
    // to explain multiplies its share by the (pessimistic) chance a base is right, so the closest contig wins most.
    private void addVotes(Map<String, SAMRecord> byContig, Map<String, Double> votesByContig)
    {
        int minDivergence = byContig.values().stream().mapToInt(ContigStatsCalculator::divergence).min().orElseThrow();

        Map<String, Double> weights = new HashMap<>();
        double weightSum = 0;
        for(Map.Entry<String, SAMRecord> entry : byContig.entrySet())
        {
            double weight = Math.pow(mCorrectBaseProbability, divergence(entry.getValue()) - minDivergence);
            weights.put(entry.getKey(), weight);
            weightSum += weight;
        }

        for(Map.Entry<String, Double> entry : weights.entrySet())
        {
            votesByContig.merge(entry.getKey(), entry.getValue() / weightSum, Double::sum);
        }
    }

    // A read aligning to >= 2 contigs is contested between strains. If one contig strictly beats the rest, the read's
    // margin (runner-up divergence minus best) is credited to it. A read that ties for best credits no contig, since
    // none stands out among indistinguishable strains.
    private static void addMargin(Map<String, SAMRecord> byContig, Map<String, List<Integer>> marginsByContig)
    {
        if(byContig.size() < 2)
        {
            return;
        }

        List<Map.Entry<String, SAMRecord>> ranked = byContig.entrySet().stream()
                .sorted(Comparator.<Map.Entry<String, SAMRecord>>comparingInt(entry -> divergence(entry.getValue())))
                .toList();

        int margin = divergence(ranked.get(1).getValue()) - divergence(ranked.get(0).getValue());
        if(margin == 0)
        {
            return;
        }

        marginsByContig.computeIfAbsent(ranked.get(0).getKey(), k -> new ArrayList<>()).add(margin);
    }

    private static ContigStats computeContig(
            String contig, int length, Collection<SAMRecord> reads, double readVotes, @Nullable List<Integer> margins)
    {
        int[] depth = new int[length];
        int[] scores = new int[reads.size()];
        int readIndex = 0;
        for(SAMRecord read : reads)
        {
            scores[readIndex++] = alignerScore(read);
            for(AlignmentBlock block : read.getAlignmentBlocks())
            {
                int start = max(0, block.getReferenceStart() - 1);   // aligned blocks are 1-based
                int end = min(length, block.getReferenceStart() - 1 + block.getLength());
                for(int i = start; i < end; ++i)
                {
                    ++depth[i];
                }
            }
        }

        int coveredBases = 0;
        for(int d : depth)
        {
            if(d > 0)
            {
                ++coveredBases;
            }
        }

        int readsBestInRivals = margins == null ? 0 : margins.size();
        Optional<SummaryStats> marginSummary = margins == null ? Optional.empty() : Optional.of(SummaryStats.from(margins));

        return new ContigStats(
                contig, length, reads.size(), coveredBases, SummaryStats.from(depth), SummaryStats.from(scores),
                readVotes, readsBestInRivals, marginSummary);
    }

    private static SAMRecord better(SAMRecord a, SAMRecord b)
    {
        int scoreA = alignerScore(a);
        int scoreB = alignerScore(b);
        if(scoreA != scoreB)
        {
            return scoreA > scoreB ? a : b;
        }
        // Deterministic tie break on location.
        return a.getAlignmentStart() <= b.getAlignmentStart() ? a : b;
    }

    private static int alignerScore(SAMRecord record)
    {
        Integer score = record.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE);
        if(score == null)
        {
            throw new IllegalStateException("aligned read missing alignment score on contig: " + record.getReferenceName());
        }
        return score;
    }

    // How poorly a contig explains the read: mismatches and indels (the NM tag) plus clipped bases. Clips count
    // because they leave read bases unexplained, so a clipped alignment must not look closer than a full-length one.
    private static int divergence(SAMRecord record)
    {
        Integer editDistance = record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE);
        if(editDistance == null)
        {
            throw new IllegalStateException("aligned read missing edit distance (NM) on contig: " + record.getReferenceName());
        }

        int clippedBases = leftClipLength(record.getCigar()) + rightClipLength(record.getCigar());
        return editDistance + clippedBases;
    }
}
