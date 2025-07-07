package com.hartwig.hmftools.bamtools.remapper;

import static com.hartwig.hmftools.bamtools.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.bamtools.common.CommonUtils.BT_LOGGER;
import static com.hartwig.hmftools.common.bamops.BamToolName.fromPath;
import static com.hartwig.hmftools.common.immune.ImmuneRegions.getHlaRegions;
import static com.hartwig.hmftools.common.perf.PerformanceCounter.runTimeMinsStr;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.bam.FastBamWriter;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class AltContigRemapper
{
    private final AltContigRemapperConfig mConfig;

    public AltContigRemapper(final AltContigRemapperConfig config)
    {
        mConfig = config;
    }

    public void run()
    {
        long startTimeMs = System.currentTimeMillis();
        BT_LOGGER.info("starting Alt-Contig remapper");

        try(SamReader samReader = SamReaderFactory.makeDefault().open(new File(mConfig.OrigBamFile)))
        {

            SAMFileHeader fileHeader = samReader.getFileHeader();
            SAMFileHeader newHeader = fileHeader.clone();

            String interimOutputFileName = mConfig.OutputFile + ".unsorted";
            File interimOutputFile = new File(interimOutputFileName);
            try(SAMFileWriter bamWriter = new FastBamWriter(newHeader, interimOutputFileName))
            {
                final BwaHlaRecordPairAligner aligner = new BwaHlaRecordPairAligner(mConfig.pairAligner(), newHeader, mConfig.RefGenVersion);
                HlaTransformer transformer = new HlaTransformer(aligner);

                if(mConfig.SliceHlaRegionsOnly)
                {
                    List<ChrBaseRegion> hlaRegions = new ArrayList<>();
                    fileHeader.getSequenceDictionary().getSequences().forEach(record ->
                    {
                        if(record.getSequenceName().toLowerCase().startsWith("hla"))
                        {
                            ChrBaseRegion baseRegion =  new ChrBaseRegion(record.getSequenceName(), record.getStart(), record.getEnd());
                            hlaRegions.add(baseRegion);
                        }
                    });

                    hlaRegions.addAll(getHlaRegions(mConfig.RefGenVersion));
                    BamSlicer bamSlicer = new BamSlicer(0, false, false, false);
                    hlaRegions.forEach(region ->
                            bamSlicer.slice(samReader, region, record -> transformer.process(record).forEach(bamWriter::addAlignment)));
                }
                else
                {
                    samReader.forEach(record -> transformer.process(record).forEach(bamWriter::addAlignment));
                }

                BT_LOGGER.info("input BAM processed, reads({}) hlaReads({})",
                        transformer.totalReadsProcessed(), transformer.hlaRecordsProcessed());

                // Deal with any unmatched reads.
                List<SAMRecord> unmatched = transformer.unmatchedRecords();
                if(!unmatched.isEmpty())
                {
                    BT_LOGGER.warn("HLA contig records unmatched({})", unmatched);
                    SingleRecordAligner unmatchedRecordsAligner = mConfig.singleRecordAligner(newHeader);
                    unmatched.forEach(samRecord -> unmatchedRecordsAligner.alignSequence(samRecord).forEach(bamWriter::addAlignment));
                }
                else
                {
                    BT_LOGGER.info("no HLA contig records were unmatched");
                }
            }
            catch(Exception e)
            {
                BT_LOGGER.error("Error processing or writing records", e);
                throw new RuntimeException(e);
            }

            BT_LOGGER.debug("BAM writer closed");

            // If the samtools path has been provided, sort the output. Else simply rename the unsorted file.
            if(mConfig.BamToolPath != null)
            {
                BT_LOGGER.debug("sorting output");
                writeSortedBam(interimOutputFileName, mConfig.OutputFile, mConfig.BamToolPath, mConfig.Threads);
                BT_LOGGER.info("sorting complete");
            }
            else
            {
                File outputFile = new File(mConfig.OutputFile);
                boolean renamed = interimOutputFile.renameTo(outputFile);
                if(!renamed)
                {
                    BT_LOGGER.warn("Could not rename " + interimOutputFile + " to " + outputFile);
                }
            }
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
        BT_LOGGER.info("remapping complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private static void writeSortedBam(final String unsortedBam, final String sortedBam, final String bamToolPath, final int threads)
    {
        if(bamToolPath == null)
        {
            return;
        }

        BT_LOGGER.info("writing sorted BAM: {}", sortedBam);

        BamToolName toolName = fromPath(bamToolPath);

        boolean success = BamOperations.sortBam(toolName, bamToolPath, unsortedBam, sortedBam, threads);

        if(success && toolName == BamToolName.SAMTOOLS)
        {
            success = BamOperations.indexBam(toolName, bamToolPath, sortedBam, threads);
        }

        if(success)
        {
            try
            {
                Files.deleteIfExists(Paths.get(unsortedBam));
            }
            catch(IOException e)
            {
                BT_LOGGER.error("error deleting interim file: {}", e.toString());
            }
        }
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        AltContigRemapperConfig.addConfig(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        AltContigRemapperConfig config = new AltContigRemapperConfig(configBuilder);
        new AltContigRemapper(config).run();
    }
}
