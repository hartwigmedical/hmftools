package com.hartwig.hmftools.virusdetect;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.function.Function.identity;
import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.stream.StreamSupport;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Per-contig support statistics over an aligned BAM: each read counted once per contig, by its single best
// alignment there.
public class ContigStatsCalculator
{
    public Map<String, ContigStats> compute(String bamFile, ViralReference reference)
    {
        return readBestAlignments(bamFile).entrySet().stream().collect(toMap(
                Map.Entry::getKey,
                entry -> computeContig(entry.getKey(), reference.contig(entry.getKey()).length(), entry.getValue().values())));
    }

    private static Map<String, Map<String, SAMRecord>> readBestAlignments(String bamFile)
    {
        try(SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile)))
        {
            return StreamSupport.stream(reader.spliterator(), false)
                    .filter(record -> !record.getReadUnmappedFlag())
                    .collect(groupingBy(SAMRecord::getReferenceName,
                            toMap(SAMRecord::getReadName, identity(), ContigStatsCalculator::better)));
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to read aligned BAM", e);
        }
    }

    private static ContigStats computeContig(String contig, int length, Collection<SAMRecord> reads)
    {
        int[] depth = new int[length];
        long scoreSum = 0;
        for(SAMRecord read : reads)
        {
            scoreSum += alignerScore(read);
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
        int minDepth = length > 0 ? Integer.MAX_VALUE : 0;
        int maxDepth = 0;
        long depthSum = 0;
        for(int d : depth)
        {
            if(d > 0)
            {
                ++coveredBases;
            }
            minDepth = min(minDepth, d);
            maxDepth = max(maxDepth, d);
            depthSum += d;
        }

        double meanDepth = length > 0 ? (double) depthSum / length : 0;
        double meanScore = (double) scoreSum / reads.size();
        return new ContigStats(contig, length, reads.size(), coveredBases, minDepth, maxDepth, meanDepth, meanScore);
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
}
