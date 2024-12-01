package com.hartwig.hmftools.redux.write;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.BAM_INDEX_EXTENSION;
import static com.hartwig.hmftools.redux.ReduxConfig.RD_LOGGER;
import static com.hartwig.hmftools.redux.ReduxConfig.humanChromosomeRegions;
import static com.hartwig.hmftools.redux.common.Constants.FILE_ID;
import static com.hartwig.hmftools.redux.common.ReadUnmapper.unmapMateAlignment;
import static com.hartwig.hmftools.redux.common.ReadUnmapper.unmapReadAlignment;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.bamops.BamOperations;
import com.hartwig.hmftools.common.bamops.BamToolName;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.region.SpecificRegions;
import com.hartwig.hmftools.redux.BamReader;
import com.hartwig.hmftools.redux.PartitionReader;
import com.hartwig.hmftools.redux.ReduxConfig;
import com.hartwig.hmftools.redux.common.FragmentStatus;
import com.hartwig.hmftools.redux.common.Statistics;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;

public class FinalBamProcessor
{
    private final ReduxConfig mConfig;
    private final FileWriterCache mFileWriterCache;
    private final BamToolName mBamToolName;

    private String mMiscBamFilename;
    private String mCombinedSortedBamFilename;
    private Statistics mStats;
    private long mUnmappedReadCount;
    private long mAltContigReadCount;

    public FinalBamProcessor(final ReduxConfig config, final FileWriterCache fileWriterCache)
    {
        mConfig = config;
        mFileWriterCache = fileWriterCache;
        mBamToolName = BamToolName.fromPath(mConfig.BamToolPath);

        mMiscBamFilename = "";
        mCombinedSortedBamFilename = "";
        mUnmappedReadCount = 0;
        mAltContigReadCount = 0;
    }

    public Statistics statistics() { return mStats; }
    public long unmappedReadCount() { return mUnmappedReadCount; }
    public long altContigReadCount() { return mAltContigReadCount; }

    public boolean run()
    {
        // first sort the unsorted BAM with the reads which have been unmapped
        if(!processUnmappedReads())
            return false;

        if(!concatenateSortedBams())
            return false;

        // final merge with concatenated partition sorted BAMs
        if(!buildFinalBam())
            return false;

        deleteInterimBams();
        return true;
    }

    private boolean processUnmappedReads()
    {
        RD_LOGGER.debug("sorting unmapped reads BAM");

        if(!mFileWriterCache.sortUnsortedBam())
            return false;

        String sortedUnmappedBamFilename = mFileWriterCache.sortedBamFilename();
            mMiscBamFilename = formInterimBamFilename("unmapped_alt_sorted");

        SAMFileWriter samFileWriter = mFileWriterCache.initialiseSamFileWriter(mMiscBamFilename, true);

        BamWriter bamWriter = new BamWriterNoSync(
                mMiscBamFilename, mConfig, mFileWriterCache.readDataWriter(), samFileWriter, null, null);

        SpecificRegions specificRegions = new SpecificRegions(); // or use the configured ones? mConfig.SpecificChrRegions
        List<ChrBaseRegion> allRegions = humanChromosomeRegions(specificRegions, mConfig.RefGenVersion);

        RD_LOGGER.debug("reprocessing unmapped reads");

        PartitionReader partitionReader = new PartitionReader(mConfig, allRegions, List.of(sortedUnmappedBamFilename), bamWriter, null);
        partitionReader.disableUnmapping();
        partitionReader.call();

        mStats = partitionReader.statistics();

        RD_LOGGER.debug("reprocessing unmapped reads");
        writeUnmappedReads(bamWriter);

        bamWriter.close();

        // index this misc unmapped and alt contig BAM
        return BamOperations.indexBam(mBamToolName, mConfig.BamToolPath, mMiscBamFilename, mConfig.Threads);
    }

    private String formInterimBamFilename(final String interimBamId)
    {
        String filename = mConfig.OutputDir + mConfig.SampleId + "." + FILE_ID;

        if(mConfig.OutputId != null)
            filename += "." + mConfig.OutputId;

        if(!interimBamId.isEmpty())
            filename += "." + interimBamId;

        filename += BAM_EXTENSION;
        return filename;
    }

    private void writeUnmappedReads(final BamWriter bamWriter)
    {
        if(mConfig.SpecificChrRegions.hasFilters() || !mConfig.WriteBam)
            return;

        BamReader bamReader = new BamReader(mConfig.BamFiles, mConfig.RefGenomeFile);

        AtomicLong unmappedCount = new AtomicLong();
        AtomicLong nonHumanContigCount = new AtomicLong();

        // do the same for non-human contigs
        bamReader.queryNonHumanContigs((final SAMRecord record) ->
        {
            processNonHumanContigReads(record, bamWriter);
            nonHumanContigCount.incrementAndGet();
        });

        bamReader.queryUnmappedReads((final SAMRecord record) ->
        {
            bamWriter.writeRead(record, FragmentStatus.UNSET);
            unmappedCount.incrementAndGet();
        });

        mUnmappedReadCount = unmappedCount.get();
        mAltContigReadCount = nonHumanContigCount.get();

        if(unmappedCount.get() > 0 || nonHumanContigCount.get() > 0)
        {
            RD_LOGGER.debug("wrote unmapped({}) altContig({}) reads", unmappedCount, nonHumanContigCount);
        }
    }

    private void processNonHumanContigReads(final SAMRecord record, final BamWriter bamWriter)
    {
        // if these have a mate in a human chromosome, then they have been unmapped in that read, so do so here as well
        if(record.getReadPairedFlag() && !record.getMateUnmappedFlag() && HumanChromosome.contains(record.getMateReferenceName()))
        {
            if(record.getSupplementaryAlignmentFlag())
                return; // drop as per standard logic

            boolean mateUnmapped = mConfig.UnmapRegions.mateInUnmapRegion(record);

            // if the human-chromosome mate was unmapped (ie in an unmap region), then this read should also now be unmapped
            // otherwise it should be unmapped but leave its mate attributes as-is
            unmapReadAlignment(record, mateUnmapped, mateUnmapped);

            if(mateUnmapped)
            {
                unmapMateAlignment(record, false, true);
            }
        }

        bamWriter.writeRead(record, FragmentStatus.UNSET);
    }

    private boolean concatenateSortedBams()
    {
        mCombinedSortedBamFilename = formInterimBamFilename("combined_sorted");

        List<String> orderPartitionBams = mFileWriterCache.bamWriters().stream()
                .filter(x -> x.isSorted()).map(x -> x.filename()).collect(Collectors.toList());

        RD_LOGGER.debug("concatenating {} sorted partition BAMs", orderPartitionBams.size());

        if(!BamOperations.concatenateBams(mBamToolName, mConfig.BamToolPath, mCombinedSortedBamFilename, orderPartitionBams, mConfig.Threads))
            return false;

        if(!BamOperations.indexBam(mBamToolName, mConfig.BamToolPath, mCombinedSortedBamFilename, mConfig.Threads))
            return false;

        return true;
    }

    private boolean buildFinalBam()
    {
        String finalBamFilename = mConfig.OutputBam != null ? mConfig.OutputBam : formInterimBamFilename("");

        RD_LOGGER.debug("writing final BAM({})", finalBamFilename);

        List<String> inputBams = List.of(mCombinedSortedBamFilename, mMiscBamFilename);

        if(!BamOperations.mergeBams(mBamToolName, mConfig.BamToolPath, finalBamFilename, inputBams, mConfig.Threads))
            return false;

        return BamOperations.indexBam(mBamToolName, mConfig.BamToolPath, finalBamFilename, mConfig.Threads);
    }

    private void deleteInterimBams()
    {
        if(mConfig.KeepInterimBams)
            return;

        List<String> interimBams = Lists.newArrayList();

        interimBams.add(mFileWriterCache.sortedBamFilename());
        interimBams.add(mMiscBamFilename);
        interimBams.add(mCombinedSortedBamFilename);

        mFileWriterCache.bamWriters().forEach(x -> interimBams.add(x.mFilename));

        try
        {
            for(String filename : interimBams)
            {
                Files.deleteIfExists(Paths.get(filename));
                Files.deleteIfExists(Paths.get(filename + BAM_INDEX_EXTENSION));
            }
        }
        catch(IOException e)
        {
            RD_LOGGER.error("error deleting interim bams: {}", e.toString());
        }
    }
}
