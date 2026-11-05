package com.hartwig.hmftools.viridian.profiling.align;

import static java.util.stream.Collectors.toMap;

import static com.hartwig.hmftools.common.bam.CigarUtils.leftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightClipLength;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.viridian.common.Utils.getRequiredSamTag;
import static com.hartwig.hmftools.viridian.common.Utils.parseMdTagMismatchPositions;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import com.hartwig.hmftools.viridian.profiling.common.ReadId;
import com.hartwig.hmftools.viridian.reference.VirusGenome;
import com.hartwig.hmftools.viridian.reference.VirusReference;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Alignments of reads to (possibly multiple) virus genomes (via BWA-MEM -a mode).
// Alignments straddling a virus genome's origin are excluded here to avoid circular genome linearization artifacts.
public record AllGenomesReadAlignments(
        ViralReadAlignmentStore store,
        AllGenomesReadAlignmentsMetrics metrics
)
{
    public static AllGenomesReadAlignments from(Stream<ViralReadAlignment> alignments)
    {
        Map<VirusGenome, Set<ReadId>> originClippedReadsByGenome = new HashMap<>();
        // BWA-MEM can place a read on one virus genome more than once. Only the best placement per virus genome is
        // kept, and the number of placements is counted. Reducing as alignments arrive avoids holding every alignment
        // in memory.
        Map<ReadGenome, BestAlignment> bestByReadGenome = new LinkedHashMap<>();

        alignments.forEach(alignment ->
        {
            if(alignment.clipsOverGenomeEnd())
            {
                originClippedReadsByGenome
                        .computeIfAbsent(alignment.genome(), k -> new HashSet<>())
                        .add(alignment.readId());
            }
            else
            {
                bestByReadGenome.merge(
                        new ReadGenome(alignment.readId(), alignment.genome()), new BestAlignment(alignment, 1),
                        BestAlignment::merge);
            }
        });

        List<ViralReadAlignment> retained = new ArrayList<>();
        Map<VirusGenome, List<Integer>> alignmentCountsByGenome = new HashMap<>();

        bestByReadGenome.forEach((readGenome, best) ->
        {
            retained.add(best.alignment());
            alignmentCountsByGenome
                    .computeIfAbsent(readGenome.genome(), k -> new ArrayList<>())
                    .add(best.alignmentCount());
        });

        Map<VirusGenome, Integer> originClippedReads = originClippedReadsByGenome.entrySet().stream()
                .collect(toMap(Map.Entry::getKey, entry -> entry.getValue().size()));

        AllGenomesReadAlignmentsMetrics metrics =
                new AllGenomesReadAlignmentsMetrics(originClippedReads, alignmentCountsByGenome);

        return new AllGenomesReadAlignments(new ViralReadAlignmentStore(retained), metrics);
    }

    public static AllGenomesReadAlignments load(String bamFile, VirusReference reference)
    {
        AlignmentDecoder decoder = new AlignmentDecoder(reference);

        try(SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bamFile));
                SAMRecordIterator records = reader.iterator())
        {
            return from(records.stream().filter(record -> !record.getReadUnmappedFlag()).map(decoder::decode));
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to read aligned BAM", e);
        }
    }

    // Helper for converting SAMRecord to ViralReadAlignment.
    static class AlignmentDecoder
    {
        private final VirusReference mReference;
        // A read aligns to many virus genomes and most alignments share a few CIGARs, so alignments share one read ID
        // and one CIGAR object instead of each holding its own copy.
        private final Map<String, ReadId> mReadIds = new HashMap<>();
        private final Map<String, Cigar> mCigars = new HashMap<>();

        AlignmentDecoder(VirusReference reference)
        {
            mReference = reference;
        }

        ViralReadAlignment decode(SAMRecord record)
        {
            ReadId readId = mReadIds.computeIfAbsent(record.getReadName(), ReadId::parse);
            Cigar cigar = mCigars.computeIfAbsent(record.getCigarString(), k -> record.getCigar());

            int editDistance = getRequiredSamTag(record, NUM_MUTATONS_ATTRIBUTE, SAMRecord::getIntegerAttribute);
            int alignerScore = getRequiredSamTag(record, ALIGNMENT_SCORE_ATTRIBUTE, SAMRecord::getIntegerAttribute);
            String mdTag = getRequiredSamTag(record, MISMATCHES_AND_DELETIONS_ATTRIBUTE, SAMRecord::getStringAttribute);
            int divergence = editDistance + leftClipLength(cigar) + rightClipLength(cigar);

            return new ViralReadAlignment(
                    readId, mReference.genome(record.getReferenceName()), record.getAlignmentStart(), cigar,
                    alignerScore, divergence, parseMdTagMismatchPositions(record.getAlignmentStart(), mdTag));
        }
    }

    private record ReadGenome(
            ReadId readId,
            VirusGenome genome
    )
    {
    }

    private record BestAlignment(
            ViralReadAlignment alignment,
            int alignmentCount
    )
    {
        BestAlignment merge(BestAlignment later)
        {
            ViralReadAlignment best = ViralReadAlignment.BEST_FIRST_TOTAL_ORDER.compare(later.alignment, alignment) < 0
                    ? later.alignment : alignment;
            return new BestAlignment(best, alignmentCount + later.alignmentCount);
        }
    }
}
