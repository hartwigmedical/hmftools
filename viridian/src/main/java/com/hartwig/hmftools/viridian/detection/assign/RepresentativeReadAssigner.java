package com.hartwig.hmftools.viridian.detection.assign;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.viridian.detection.align.AlignedRead;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignment;
import com.hartwig.hmftools.viridian.detection.align.ViralReadAlignments;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

// Selects each read's best alignment to the selected representative virus contig.
// Writes those alignments to a coordinate sorted, indexed BAM.
// Reads with no alignments to the representative are dropped.
public class RepresentativeReadAssigner
{
    private static final Logger LOGGER = LogManager.getLogger(RepresentativeReadAssigner.class);

    public static void assign(
            ViralReadAlignments viralAlignments, Set<ViralContig> representatives,
            String allAlignmentsBamFile, String outputBamFile)
    {
        Map<String, ViralReadAlignment> assignmentsByRead = assignReads(viralAlignments.reads(), representatives);

        LOGGER.debug(
                "Assigned {} reads, {} fit no representative",
                assignmentsByRead.size(), viralAlignments.reads().size() - assignmentsByRead.size());

        writeAssignedAlignments(allAlignmentsBamFile, outputBamFile, assignmentsByRead);
    }

    // Return value map is: read name -> selected alignment
    // Absent if the read has no alignment.
    static Map<String, ViralReadAlignment> assignReads(List<AlignedRead> reads, Set<ViralContig> representatives)
    {
        // Tie-break by contig for determinism.
        // But not expecting that a read will support multiple representative contigs.
        Comparator<ViralReadAlignment> alignmentComparator =
                ViralReadAlignment.BEST_FIT_FIRST.thenComparing(ViralReadAlignment::contig);

        Map<String, ViralReadAlignment> assignmentsByRead = new LinkedHashMap<>();
        int multiRepReads = 0;
        int multiRepTiedReads = 0;
        for(AlignedRead read : reads)
        {
            List<ViralReadAlignment> candidates = read.hits().entrySet().stream()
                    .filter(hit -> representatives.contains(hit.getKey()))
                    .map(hit -> hit.getValue().best())
                    .sorted(alignmentComparator)
                    .toList();
            if(!candidates.isEmpty())
            {
                assignmentsByRead.put(read.readName(), candidates.get(0));
            }
            if(candidates.size() > 1)
            {
                ++multiRepReads;
                if(ViralReadAlignment.BEST_FIT_FIRST.compare(candidates.get(0), candidates.get(1)) == 0)
                {
                    // Not only does this read support two representatives, but they are indistinguishable.
                    ++multiRepTiedReads;
                }
            }
        }

        // Log the cases of 1 read supporting multiple representative contigs.
        // Not expecting it, but it could occur - we need more data.
        if(multiRepReads > 0)
        {
            LOGGER.debug("{} reads fit several representatives, {} of those indistinguishably", multiRepReads, multiRepTiedReads);
        }

        return assignmentsByRead;
    }

    // Slice the BAM on the selected alignments.
    private static void writeAssignedAlignments(
            String sourceBamFile, String outputBamFile, Map<String, ViralReadAlignment> assignmentsByRead)
    {
        int writtenRecords = 0;

        try(SamReader reader = SamReaderFactory.makeDefault()
                .validationStringency(ValidationStringency.SILENT).open(new File(sourceBamFile)))
        {
            SAMFileHeader header = reader.getFileHeader().clone();
            header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

            try(SAMFileWriter writer = new SAMFileWriterFactory()
                    .setCreateIndex(true)
                    .makeBAMWriter(header, false, new File(outputBamFile)))
            {
                for(SAMRecord record : reader)
                {
                    ViralReadAlignment assigned = assignmentsByRead.get(record.getReadName());
                    if(assigned != null && isSameAlignment(assigned, record))
                    {
                        record.setHeaderStrict(header);
                        // The read keeps only this alignment here, so a secondary or supplementary mark would point at a
                        // primary not in this BAM, and viewers hide such records.
                        record.setSecondaryAlignment(false);
                        record.setSupplementaryAlignmentFlag(false);
                        writer.addAlignment(record);
                        ++writtenRecords;
                    }
                }
            }

            fixBamIndexName(outputBamFile);
        }
        catch(IOException e)
        {
            throw new RuntimeException("Failed to write representative alignment BAM", e);
        }

        // Double-check that the correct alignments were written; otherwise the downstream stats will be wrong.
        if(writtenRecords != assignmentsByRead.size())
        {
            throw new IllegalStateException(String.format(
                    "Wrote %d records for %d assigned reads", writtenRecords, assignmentsByRead.size()));
        }
    }

    private static boolean isSameAlignment(ViralReadAlignment assigned, SAMRecord record)
    {
        return assigned.alignmentStart() == record.getAlignmentStart()
                && assigned.contig().name().equals(record.getReferenceName())
                // Note CIGAR is needed too because occasionally there are multiple plausible alignments at the same ref position.
                && assigned.cigar().equals(record.getCigarString());
    }

    // htsjdk names the index after the BAM's file stem, but tooling looks for it appended to the full name.
    // So "HG002.viridian.representative.bai" is renamed to "HG002.viridian.representative.bam.bai".
    private static void fixBamIndexName(String bamFile) throws IOException
    {
        Path written = Path.of(bamFile.substring(0, bamFile.length() - BAM_EXTENSION.length()) + BAM_INDEX_EXTENSION);
        Files.move(written, Path.of(bamFile + BAM_INDEX_EXTENSION), StandardCopyOption.REPLACE_EXISTING);
    }
}
