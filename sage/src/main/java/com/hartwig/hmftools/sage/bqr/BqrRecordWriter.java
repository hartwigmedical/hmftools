package com.hartwig.hmftools.sage.bqr;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.StringJoiner;

import com.hartwig.hmftools.sage.SageConfig;

import htsjdk.samtools.SAMRecord;

public class BqrRecordWriter
{
    private final SageConfig mConfig;
    private final BufferedWriter mWriter;

    public BqrRecordWriter(final SageConfig config, final String sampleId)
    {
        mConfig = config;
        mWriter = initialiseWriter(sampleId);
    }

    public boolean enabled() { return mWriter != null; }

    public synchronized void writeRecordData(
            final SAMRecord record,
            final int position, final int readIndex, final byte ref, final byte alt, final byte[] trinucleotideContext, final byte quality)
    {
        if(mWriter == null)
            return;

        /*
        input: trinucleotide, alt
        output: per read matching these criteria, output a row containing
        (chrom, pos, read ID, qual,
        fwd/bwd strand, fragment orientation,
        skip aggregate: number of reads in fragment per strand,
        distance from 3' end of read,
        distance from 5' end of read,
        skip aggregate: number of errors at this site on fwd strand with any qual,
        skip aggregate: number of errors at this site on bwd strand with any qual)
        */

        try
        {
            mWriter.write(format("%c\t%c\t%s", (char)ref, (char)alt, new String(trinucleotideContext)));
            mWriter.write(format("\t%s\t%d\t%d\t%d", record.getReferenceName(), position, readIndex, quality));

            boolean isNegStrand = record.getReadNegativeStrandFlag();
            int readLength = record.getReadBases().length;
            int distanceFromFivePrimeEnd = isNegStrand ? readLength - readIndex - 1 : readIndex;
            int distanceFromThreePrimeEnd = readLength - distanceFromFivePrimeEnd;

            mWriter.write(format("\t%c\t%s\t%d\t%d",
                    orientationChar(record.getReadNegativeStrandFlag()), fragmentOrientationStr(record),
                    distanceFromFivePrimeEnd, distanceFromThreePrimeEnd));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write BQR read data", e.toString());
        }
    }

    private static String fragmentOrientationStr(final SAMRecord record)
    {
        if(record.getFirstOfPairFlag())
            return format("%c1%c2", orientationChar(record.getReadNegativeStrandFlag()), orientationChar(record.getMateNegativeStrandFlag()));
        else
            return format("%c1%c2", orientationChar(record.getMateNegativeStrandFlag()), orientationChar(record.getReadNegativeStrandFlag()));
    }

    private static char orientationChar(boolean negStrand) { return negStrand ? 'R' : 'F'; }

    public void close() { closeBufferedWriter(mWriter); }

    private BufferedWriter initialiseWriter(final String sampleId)
    {
        if(!mConfig.QualityRecalibration.WriteReads)
            return null;

        try
        {
            String outputFile = mConfig.outputDir() + sampleId + ".bqr_reads.tsv";;
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add("Ref").add("Alt").add("TriNuc");
            sj.add("Chromosome").add("Position").add("ReadIndex").add("Quality");
            sj.add("ReadOrient").add("FragOrient").add("FivePrimeDist").add("ThreePrimeDist");

            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            e.printStackTrace();
            SG_LOGGER.error("failed to write BQR read data", e.toString());
            System.exit(1);
        }

        return null;
    }
}
