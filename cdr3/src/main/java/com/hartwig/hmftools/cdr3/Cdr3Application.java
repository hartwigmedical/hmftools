package com.hartwig.hmftools.cdr3;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.checkCreateOutputDir;

import java.io.File;
import java.io.IOException;
import java.time.Duration;
import java.time.Instant;
import java.time.LocalDateTime;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.function.BiConsumer;
import java.util.stream.Collectors;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParametersDelegate;
import com.beust.jcommander.UnixStyleUsageFormatter;
import com.hartwig.hmftools.cdr3.layout.ReadLayout;
import com.hartwig.hmftools.cdr3.layout.VDJCandidate;
import com.hartwig.hmftools.common.codon.Codons;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.utils.config.DeclaredOrderParameterComparator;
import com.hartwig.hmftools.common.utils.config.LoggingOptions;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.ReferenceSource;
import kotlin.ranges.IntRange;

public class Cdr3Application
{
    public static final Logger sLogger = LogManager.getLogger(Cdr3Application.class);
    
    // add the options
    @ParametersDelegate
    private final Cdr3Params mParams = new Cdr3Params();

    // add to the logging options
    @ParametersDelegate
    private final LoggingOptions mLoggingOptions = new LoggingOptions();

    public int run() throws IOException, InterruptedException
    {
        mLoggingOptions.setLogLevel();

        VersionInfo versionInfo = new VersionInfo("cdr3.version");
        sLogger.info("CDR3 version: {}, build time: {}", versionInfo.version(), versionInfo.buildTime());

        if(!mParams.isValid())
        {
            sLogger.error(" invalid config, exiting");
            return 1;
        }

        checkCreateOutputDir(mParams.OutputDir);

        Instant start = Instant.now();

        VJGeneStore vjGeneStore = new Cdr3GeneLoader(mParams.refGenomeVersion);
        Cdr3ReadScreener readProcessor = new Cdr3ReadScreener(vjGeneStore, mParams.MaxAnchorAlignDistance, CiderConstants.MIN_CANDIDATE_READ_ANCHOR_OVERLAP);
        readBamFile(readProcessor, vjGeneStore);

        //VJReadJoiner vjReadJoiner = new VJReadJoiner(readProcessor.getCdr3ReadMatches());

        Instant finish = Instant.now();
        long seconds = Duration.between(start, finish).getSeconds();
        sLogger.info("CDR3 run complete, time taken: {}m {}s", seconds / 60, seconds % 60);

        return 0;
    }

    public void readBamFile(Cdr3ReadScreener readProcessor, VJGeneStore vjGeneStore) throws InterruptedException, IOException
    {
        final SamReaderFactory readerFactory = readerFactory(mParams);

        // create a thread
        final BlockingQueue<SAMRecord> samRecordQueue = new LinkedBlockingQueue<>();

        final Thread readProcessorThread = new Thread(() ->
        {
            while (true)
            {
                try
                {
                    SAMRecord record = samRecordQueue.take();

                    if (record.getHeader() == null)
                        // indicates finish
                        return;
                    readProcessor.processSamRecord(record);
                } catch (InterruptedException e)
                {
                    break;
                }
            }
        });

        readProcessorThread.start();

        BiConsumer<GenomeRegion, SAMRecord> lociBamRecordHander = (GenomeRegion genomeRegion, SAMRecord samRecord) ->
        {
            samRecordQueue.add(samRecord);
        };

        Collection<GenomeRegion> genomeRegions = vjGeneStore.getVJAnchorReferenceLocations().stream()
                .map(o -> GenomeRegions.create(o.getChromosome(), o.getStart() - mParams.MaxAnchorAlignDistance, o.getEnd() + mParams.MaxAnchorAlignDistance))
                .collect(Collectors.toList());

        AsyncBamReader.processBam(mParams.BamPath, readerFactory, genomeRegions, lociBamRecordHander, mParams.ThreadCount, 0);

        // a bit hacky to tell consumer this is the end
        samRecordQueue.add(new SAMRecord(null));

        readProcessorThread.join();

        sLogger.info("found {} VJ read records", readProcessor.getAllMatchedReads().size());

        //String readTsvFile = Cdr3ReadTsvWriter.generateFilename(mParams.OutputDir, mParams.SampleId);
        //Cdr3ReadTsvWriter.write(readTsvFile, readProcessor.getVJReadCandidates());

        // now build the consensus overlay sequences
        //var geneTypes = new VJGeneType[] { VJGeneType.IGHV, VJGeneType.IGHJ };
        var geneTypes = VJGeneType.values();
        Map<VJGeneType, List<ReadLayout>> layoutMap = new HashMap<>();

        for (VJGeneType geneType : geneTypes)
        {
            List<VJReadCandidate> readCandidates = readProcessor.getVJReadCandidates().values().stream()
                    .filter(o -> o.getVjGeneType() == geneType).collect(Collectors.toList());

            for (var read : readCandidates)
            {
                if ((read.getAnchorOffsetEnd() - read.getAnchorOffsetStart()) > 30)
                {
                    sLogger.info("anchor length > 30: read: {}, cigar: {}", read.getRead(), read.getRead().getCigarString());
                }
            }

            List<ReadLayout> readLayouts = VJReadLayoutAdaptor.buildOverlays(geneType, readCandidates,
                    mParams.MinBaseQuality, 20, 1.0, mParams.numBasesToTrim);

            // now log the sequences
            for (ReadLayout layout : readLayouts)
            {
                List<VJReadCandidate> overlayReads =
                        layout.getReads().stream().map(VJReadLayoutAdaptor::toReadCandidate).collect(Collectors.toList());
                int numSplitReads = (int) overlayReads.stream().filter(o -> Math.max(o.getLeftSoftClip(), o.getRightSoftClip()) > 5).count();
                int anchorLength =
                        overlayReads.stream().map(o -> o.getAnchorOffsetEnd() - o.getAnchorOffsetStart()).max(Integer::compareTo).orElse(0);
                sLogger.info("Overlay type: {}， read count: {}, split read count: {}, anchor length: {}",
                        geneType, layout.getReads().size(), numSplitReads, anchorLength);

                // get the sequence, remember aligned position is the anchor start
                String anchor = VJReadLayoutAdaptor.getAnchorSequence(geneType, layout);
                String cdr3 = VJReadLayoutAdaptor.getCdr3Sequence(geneType, layout);
                String anchorSupport = VJReadLayoutAdaptor.getAnchorSupport(geneType, layout);
                String cdr3Support = VJReadLayoutAdaptor.getCdr3Support(geneType, layout);

                if (geneType.getVj() == VJ.V)
                {
                    sLogger.info("V sequence: {}-{}", anchor, cdr3);
                    sLogger.info("V support:  {}-{}", anchorSupport, cdr3Support);
                    sLogger.info("V AA seq:  {}-{}", Codons.aminoAcidFromBases(anchor), Codons.aminoAcidFromBases(cdr3));
                } else if (geneType.getVj() == VJ.J)
                {
                    sLogger.info("J sequence: {}-{}", cdr3, anchor);
                    sLogger.info("J support:  {}-{}", cdr3Support, anchorSupport);
                    sLogger.info("J AA seq:  {}-{}", Codons.aminoAcidFromBases(cdr3), Codons.aminoAcidFromBases(anchor));
                }
            }

            layoutMap.put(geneType, readLayouts);
        }

        VJReadLayoutFile.writeLayouts(mParams.OutputDir, mParams.SampleId, layoutMap, 5);

        // now try to find the VDJs
        VJLayoutAligner vjLayoutAligner = new VJLayoutAligner(20, 1.0);
        List<VDJCandidate> vdjCandidates = vjLayoutAligner.findAlignedVJs(layoutMap.get(VJGeneType.IGHV), layoutMap.get(VJGeneType.IGHJ));

        // now print them out
        for (VDJCandidate vdj : vdjCandidates)
        {
            // we want to use the indices to work where things are
            String vdjSeq = vdj.getVdjSequence();
            IntRange vAnchorRange = VJReadLayoutAdaptor.getAnchorRange(VJGeneType.IGHV, vdj.getVLayout());

            // we need to shift it
            IntRange jAnchorRange = VJReadLayoutAdaptor.getAnchorRange(VJGeneType.IGHJ, vdj.getJLayout());
            int jSeqOffset = vdjSeq.length() - vdj.getJLayout().consensusSequence().length();
            jAnchorRange = new IntRange(jAnchorRange.getStart() + jSeqOffset, jAnchorRange.getEndInclusive() + jSeqOffset);

            String vAnchor = vdjSeq.substring(vAnchorRange.getStart(), vAnchorRange.getEndInclusive() + 1);
            String cdr = vdjSeq.substring(vAnchorRange.getEndInclusive() + 1, jAnchorRange.getStart());
            String jAnchor = vdjSeq.substring(jAnchorRange.getStart(), jAnchorRange.getEndInclusive() + 1);

            // now split it up into VDJ
            sLogger.debug("v and j overlap");
            sLogger.debug("v: {}", vdj.getVLayout().consensusSequence());
            sLogger.debug("j: {}", vdj.getJLayout().consensusSequence());
            sLogger.debug("overlap: {}", vdj.getOverlapSeq());
            sLogger.debug("vdj: {}-{}-{}", vAnchor, cdr, jAnchor);
            sLogger.debug("vdj AA: {}-{}-{}", Codons.aminoAcidFromBases(vAnchor), Codons.aminoAcidFromBases(cdr), Codons.aminoAcidFromBases(jAnchor));
        }

        String vdjTsvFile = Cdr3SequenceTsvWriter.generateFilename(mParams.OutputDir, mParams.SampleId);
        Cdr3SequenceTsvWriter.write(vdjTsvFile, vdjCandidates);

        List<VDJSequence> vdjSequences = VDJSequenceBuilder.buildVDJSequences(layoutMap, vjGeneStore, (byte)mParams.MinBaseQuality);

        VDJSequenceTsvWriter.writeVDJSequences(mParams.OutputDir, mParams.SampleId, vdjSequences);

        writeCdr3Bam(readProcessor.getAllMatchedReads());
    }

    public void writeCdr3Bam(Collection<SAMRecord> samRecords) throws IOException
    {
        if (!mParams.writeFilteredBam)
            return;

        String outBamPath = mParams.OutputDir + "/" + mParams.SampleId + ".cider.bam";

        final SamReaderFactory readerFactory = readerFactory(mParams);
        SAMFileHeader samFileHeader;
        try (SamReader samReader = readerFactory.open(new File(mParams.BamPath)))
        {
            samFileHeader = samReader.getFileHeader();
        }
        try (SAMFileWriter bamFileWriter = (new SAMFileWriterFactory()).makeBAMWriter(
                samFileHeader, false, new File(outBamPath)))
        {
            for (SAMRecord r : samRecords)
            {
                bamFileWriter.addAlignment(r);
            }
        }
    }

    private static SamReaderFactory readerFactory(final Cdr3Params params)
    {
        final SamReaderFactory readerFactory = SamReaderFactory.make().validationStringency(params.Stringency);
        if(params.RefGenomePath != null)
        {
            return readerFactory.referenceSource(new ReferenceSource(new File(params.RefGenomePath)));
        }
        return readerFactory;
    }

    public static void main(final String... args) throws IOException, InterruptedException
    {
        sLogger.info("{}", LocalDateTime.now());
        sLogger.info("args: {}", String.join(" ", args));

        Cdr3Application cdr3Application = new Cdr3Application();
        JCommander commander = JCommander.newBuilder()
                .addObject(cdr3Application)
                .build();

        // use unix style formatter
        commander.setUsageFormatter(new UnixStyleUsageFormatter(commander));
        // help message show in order parameters are declared
        commander.setParameterDescriptionComparator(new DeclaredOrderParameterComparator(Cdr3Application.class));

        try
        {
            commander.parse(args);
        }
        catch (com.beust.jcommander.ParameterException e)
        {
            System.out.println("Unable to parse args: " + e.getMessage());
            commander.usage();
            System.exit(1);
        }

        System.exit(cdr3Application.run());
    }
}
