package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.Nullable;

// Scans the tumor BAM/CRAM once and writes candidate viral reads single-end to a FASTA. Duplicate
// and secondary/supplementary alignments are skipped so each read is emitted once; the mate number
// is appended to the FASTA id to keep the two reads of a pair distinct.
public class CandidateReadExtractor
{
    @Nullable
    private final String mRefGenomeFile; // required only for CRAM decode
    private final CandidateReadFilter mFilter;

    private static final Logger LOGGER = LogManager.getLogger(CandidateReadExtractor.class);

    public CandidateReadExtractor(@Nullable String refGenomeFile, CandidateReadFilter filter)
    {
        mRefGenomeFile = refGenomeFile;
        mFilter = filter;
    }

    public int extractToFasta(String tumorBam, String outputFasta)
    {
        SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
        if(mRefGenomeFile != null)
        {
            factory = factory.referenceSequence(new File(mRefGenomeFile));
        }

        int candidateCount = 0;
        try(SamReader reader = factory.open(new File(tumorBam));
                BufferedWriter writer = createBufferedWriter(outputFasta))
        {
            for(SAMRecord record : reader)
            {
                if(record.getDuplicateReadFlag() || record.isSecondaryOrSupplementary())
                {
                    continue;
                }

                if(mFilter.isCandidate(record))
                {
                    writeFasta(writer, record);
                    ++candidateCount;
                }
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException("failed to extract candidate reads", e);
        }

        LOGGER.info("extracted {} candidate reads to {}", candidateCount, outputFasta);
        return candidateCount;
    }

    private static void writeFasta(BufferedWriter writer, SAMRecord record) throws IOException
    {
        writer.write(">" + fastaId(record));
        writer.newLine();
        writer.write(record.getReadString());
        writer.newLine();
    }

    private static String fastaId(SAMRecord record)
    {
        if(record.getReadPairedFlag())
        {
            return record.getReadName() + (record.getFirstOfPairFlag() ? "/1" : "/2");
        }
        return record.getReadName();
    }
}
