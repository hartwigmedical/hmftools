package com.hartwig.hmftools.svprep.tools;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConfig.SAMPLE;

import java.io.File;
import java.io.IOException;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class BamSvSlicer
{
    private final String mVcfFilename;
    private final String mBamFile;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;
    private SAMFileWriter mWriter;

    private int mBamRecords;

    private static final String VCF_FILE = "vcf_file";
    private static final String BAM_FILE = "bam_file";

    public BamSvSlicer(final CommandLine cmd)
    {
        mVcfFilename = cmd.getOptionValue(VCF_FILE);
        mBamFile = cmd.getOptionValue(BAM_FILE);

        String sampleId = cmd.getOptionValue(SAMPLE);
        String outputDir = parseOutputDir(cmd);
        String outputId = cmd.getOptionValue(OUTPUT_ID);

        String refGenomeFile = cmd.getOptionValue(REF_GENOME);

        String outputBam = outputDir + sampleId;
        if(outputId != null)
            outputBam += "." + outputId;

        outputBam += ".vcf_sliced.bam";

        mWriter = initialise(refGenomeFile, outputBam);
        mBamRecords = 0;

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(new File(refGenomeFile)).open(new File(mBamFile));
        mBamSlicer = new BamSlicer(0, true, true, true);
    }

    private static int MIN_DISTANCE = 3000;
    private static int POS_BUFFER = 1000;

    public void run()
    {
        if(mVcfFilename == null)
        {
            SV_LOGGER.error("missing VCF");
            System.exit(1);
        }

        final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                mVcfFilename, new VCFCodec(), false);

        // VCFHeader vcfHeader = (VCFHeader)reader.getHeader();

        int vcfCount = 0;
        int sliceCount = 0;

        try
        {
            String currentChromosome = "";
            int lastPosition = -1;
            int lastStartPosition = -1;
            int lastBamCount = 0;

            for(VariantContext variantContext : reader.iterator())
            {
                ++vcfCount;
                String chromosome = variantContext.getContig();
                int position = variantContext.getStart();

                if(!chromosome.equals(currentChromosome))
                {
                    if(lastPosition > 0)
                    {
                        sliceRegion(currentChromosome, lastStartPosition, lastPosition);
                        SV_LOGGER.debug("chromosome({}) SVs sliced, BAM records({})", chromosome, mBamRecords - lastBamCount);
                    }

                    SV_LOGGER.debug("processing chromosome({})", chromosome);
                    lastPosition = -1;
                    lastStartPosition = -1;
                    currentChromosome = chromosome;
                    lastBamCount = mBamRecords;
                }

                if(lastPosition > 0)
                {
                    int posDiff = position - lastPosition;

                    if(posDiff <= MIN_DISTANCE)
                    {
                        lastPosition = position;
                        continue;
                    }

                    sliceRegion(currentChromosome, lastStartPosition, lastPosition);
                }

                lastStartPosition = lastPosition = position;
            }

            if(lastPosition > 0)
                sliceRegion(currentChromosome, lastStartPosition, lastPosition);
        }
        catch(IOException e)
        {
            SV_LOGGER.error("error reading vcf({}): {}", mVcfFilename, e.toString());
            System.exit(1);
        }

        SV_LOGGER.info("BAM SV slice complete: SVs({}) slice regions({}) BAM records({})", vcfCount, sliceCount, mBamRecords);

        mWriter.close();
    }

    private SAMFileWriter initialise(final String refGenomeFile, final String outputFile)
    {
        SamReader samReader = SamReaderFactory.makeDefault().referenceSequence(new File(refGenomeFile)).open(new File(mBamFile));
        return new SAMFileWriterFactory().makeBAMWriter(samReader.getFileHeader(), false, new File(outputFile));
    }

    private void sliceRegion(final String chromosome, int posStart, int posEnd)
    {
        ChrBaseRegion region = new ChrBaseRegion(chromosome, posStart - POS_BUFFER, posEnd + POS_BUFFER);
        SV_LOGGER.debug("slicing region({}) length({})", region, region.baseLength());
        mBamSlicer.slice(mSamReader, Lists.newArrayList(region), this::processSamRecord);
    }

    private void processSamRecord(final SAMRecord record)
    {
        mWriter.addAlignment(record);
        ++mBamRecords;
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the sample");
        options.addOption(VCF_FILE, true, "VCF File");
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(BAM_FILE, true, "BAM file to slice");

        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        BamSvSlicer bamSvSlicer = new BamSvSlicer(cmd);
        bamSvSlicer.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
