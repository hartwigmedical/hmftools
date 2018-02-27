package com.hartwig.hmftools.bamslicer;

import static htsjdk.samtools.util.BlockCompressedFilePointerUtil.MAX_BLOCK_ADDRESS;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.slicing.Slicer;
import com.hartwig.hmftools.common.slicing.SlicerFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.BAMFileReader;
import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.DiskBasedBAMFileIndex;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.BlockCompressedFilePointerUtil;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import okhttp3.OkHttpClient;

public class BamSlicerApplication {
    private static final Logger LOGGER = LogManager.getLogger(BamSlicerApplication.class);
    private static final String SBP_ENDPOINT_URL = System.getenv("SBP_ENDPOINT_URL");
    private static final String SBP_PROFILE = "download";
    private static final int EXPIRATION_HOURS = 2;

    private static final String INPUT_MODE_S3 = "s3";
    private static final String INPUT_MODE_URL = "url";
    private static final String INPUT_MODE_FILE = "file";
    private static final String INPUT = "input";
    private static final String BUCKET = "bucket";
    private static final String INDEX = "index";
    private static final String OUTPUT = "output";
    private static final String PROXIMITY = "proximity";
    private static final String VCF = "vcf";
    private static final String BED = "bed";
    private static final String MAX_CHUNKS_IN_MEMORY = "max_chunks";
    private static final String MAX_CONCURRENT_REQUESTS = "max_concurrent_requests";

    public static void main(final String... args) throws ParseException, IOException {
        final CommandLine cmd = createCommandLine(args);
        assert cmd != null;
        //MIVO: disable default samtools buffering
        System.setProperty("samjdk.buffer_size", "0");
        if (cmd.hasOption(INPUT_MODE_FILE)) {
            sliceFromVCF(cmd);
        }
        if (cmd.hasOption(INPUT_MODE_S3)) {
            final Pair<URL, URL> urls = generateURLs(cmd);
            sliceFromURLs(urls.getKey(), urls.getValue(), cmd);
        }
        if (cmd.hasOption(INPUT_MODE_URL)) {
            final URL bamURL = new URL(cmd.getOptionValue(INPUT));
            final URL indexURL = new URL(cmd.getOptionValue(INDEX));
            sliceFromURLs(indexURL, bamURL, cmd);
        }
        LOGGER.info("Done.");
    }

    private static void sliceFromVCF(@NotNull final CommandLine cmd) throws IOException {
        final String inputPath = cmd.getOptionValue(INPUT);
        final String outputPath = cmd.getOptionValue(OUTPUT);
        final String vcfPath = cmd.getOptionValue(VCF);
        final int proximity = Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500"));
        final SamReader reader = SamReaderFactory.makeDefault().open(new File(inputPath));
        final QueryInterval[] intervals = getIntervalsFromVCF(vcfPath, reader.getFileHeader(), proximity);
        final CloseableIterator<SAMRecord> iterator = getIterator(reader, intervals);
        writeToSlice(outputPath, reader.getFileHeader(), iterator);
        reader.close();
    }

    @NotNull
    private static QueryInterval[] getIntervalsFromVCF(@NotNull final String vcfPath, @NotNull final SAMFileHeader header,
            final int proximity) {
        final File vcfFile = new File(vcfPath);
        final VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);
        final List<QueryInterval> queryIntervals = Lists.newArrayList();
        for (VariantContext variant : vcfReader) {

            queryIntervals.add(new QueryInterval(header.getSequenceIndex(variant.getContig()),
                    Math.max(0, variant.getStart() - proximity),
                    variant.getStart() + proximity));

            if (variant.getStructuralVariantType() == StructuralVariantType.BND) {
                final String call = variant.getAlternateAllele(0).getDisplayString();
                final String[] leftSplit = call.split("]");
                final String[] rightSplit = call.split("\\[");

                final String contig;
                final int position;
                if (leftSplit.length >= 2) {
                    final String[] location = leftSplit[1].split(":");
                    contig = location[0];
                    position = Integer.parseInt(location[1]);
                } else if (rightSplit.length >= 2) {
                    final String[] location = rightSplit[1].split(":");
                    contig = location[0];
                    position = Integer.parseInt(location[1]);
                } else {
                    LOGGER.error(variant.getID() + " : could not parse breakpoint");
                    continue;
                }
                queryIntervals.add(new QueryInterval(header.getSequenceIndex(contig),
                        Math.max(0, position - proximity),
                        position + proximity));
            } else {
                queryIntervals.add(new QueryInterval(header.getSequenceIndex(variant.getContig()),
                        Math.max(0, variant.getEnd() - proximity),
                        variant.getEnd() + proximity));
            }
        }

        return QueryInterval.optimizeIntervals(queryIntervals.toArray(new QueryInterval[queryIntervals.size()]));
    }

    @NotNull
    private static Pair<URL, URL> generateURLs(@NotNull final CommandLine cmd) {
        try {
            LOGGER.info("Attempting to generate S3 URLs for endpoint: {} using profile: {}", SBP_ENDPOINT_URL, SBP_PROFILE);
            final S3UrlGenerator urlGenerator = ImmutableS3UrlGenerator.of(SBP_ENDPOINT_URL, SBP_PROFILE);
            final URL indexUrl = urlGenerator.generateUrl(cmd.getOptionValue(BUCKET), cmd.getOptionValue(INDEX), EXPIRATION_HOURS);
            final URL bamUrl = urlGenerator.generateUrl(cmd.getOptionValue(BUCKET), cmd.getOptionValue(INPUT), EXPIRATION_HOURS);
            return Pair.of(indexUrl, bamUrl);
        } catch (Exception e) {
            LOGGER.error("Could not create S3 URLs. Error: {}", e.toString());
            LOGGER.error("You must run this with the sbp user or set up aws credentials and the SBP_ENDPOINT_URL environment variable");
            System.exit(1);
        }
        return null;
    }

    private static void sliceFromURLs(@NotNull final URL indexUrl, @NotNull final URL bamUrl, @NotNull final CommandLine cmd)
            throws IOException {
        final OkHttpClient httpClient = SlicerHttpClient.create(Integer.parseInt(cmd.getOptionValue(MAX_CONCURRENT_REQUESTS)));
        final String outputPath = cmd.getOptionValue(OUTPUT);
        final String bedPath = cmd.getOptionValue(BED);
        final int maxBufferSize = readMaxBufferSize(cmd);
        final File indexFile = downloadIndex(indexUrl);
        final SamReader reader = SamReaderFactory.makeDefault().open(SamInputResource.of(bamUrl).index(indexFile));
        LOGGER.info("Generating query intervals from BED file: {}", bedPath);
        final QueryInterval[] intervals = getIntervalsFromBED(bedPath, reader.getFileHeader());
        final BAMFileSpan span = bamSpanForIntervals(indexFile, reader.getFileHeader(), intervals);
        final List<Chunk> expandedChunks = expandChunks(span.getChunks());
        LOGGER.info("Generated {} query intervals which map to {} bam chunks", intervals.length, expandedChunks.size());
        final SamInputResource bamResource =
                SamInputResource.of(new CachingSeekableHTTPStream(httpClient, bamUrl, expandedChunks, maxBufferSize)).index(indexFile);
        final SamReader cachingReader = SamReaderFactory.makeDefault().open(bamResource);

        LOGGER.info("Slicing bam...");
        final CloseableIterator<SAMRecord> iterator = getIterator(cachingReader, intervals, span.toCoordinateArray());
        writeToSlice(outputPath, cachingReader.getFileHeader(), iterator);
        cachingReader.close();
        reader.close();
        indexFile.deleteOnExit();
    }

    @NotNull
    private static List<Chunk> expandChunks(@NotNull final List<Chunk> chunks) {
        final List<Chunk> result = Lists.newArrayList();
        //MIVO: add chunk for header
        final long headerEndVirtualPointer = ((long) BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE) << 16;
        result.add(new Chunk(0, headerEndVirtualPointer));
        for (final Chunk chunk : chunks) {
            final long chunkEndBlockAddress = BlockCompressedFilePointerUtil.getBlockAddress(chunk.getChunkEnd());
            final long extendedEndBlockAddress = chunkEndBlockAddress + BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE;
            final long newChunkEnd = extendedEndBlockAddress > MAX_BLOCK_ADDRESS ? MAX_BLOCK_ADDRESS : extendedEndBlockAddress;
            final long chunkEndVirtualPointer = newChunkEnd << 16;
            result.add(new Chunk(chunk.getChunkStart(), chunkEndVirtualPointer));
        }
        return Chunk.optimizeChunkList(result, 0);
    }

    @NotNull
    private static BAMFileSpan bamSpanForIntervals(@NotNull final File index, @NotNull final SAMFileHeader header,
            @NotNull final QueryInterval[] intervals) {
        final BAMIndex bamIndex = new DiskBasedBAMFileIndex(index, header.getSequenceDictionary(), false);
        return BAMFileReader.getFileSpan(intervals, bamIndex);
    }

    @NotNull
    private static File downloadIndex(@NotNull final URL indexUrl) throws IOException {
        LOGGER.info("Downloading index...");
        final ReadableByteChannel indexChannel = Channels.newChannel(indexUrl.openStream());
        final File index = File.createTempFile("tmp", ".bai");
        final FileOutputStream indexOutputStream = new FileOutputStream(index);
        indexOutputStream.getChannel().transferFrom(indexChannel, 0, Long.MAX_VALUE);
        indexOutputStream.close();
        indexChannel.close();
        return index;
    }

    @NotNull
    private static QueryInterval[] getIntervalsFromBED(@NotNull final String bedPath, @NotNull final SAMFileHeader header)
            throws IOException {
        final Slicer bedSlicer = SlicerFactory.fromBedFile(bedPath);
        final List<QueryInterval> queryIntervals = Lists.newArrayList();
        for (final GenomeRegion region : bedSlicer.regions()) {
            queryIntervals.add(new QueryInterval(header.getSequenceIndex(region.chromosome()), (int) region.start(), (int) region.end()));
        }
        return QueryInterval.optimizeIntervals(queryIntervals.toArray(new QueryInterval[queryIntervals.size()]));
    }

    @NotNull
    private static CloseableIterator<SAMRecord> getIterator(@NotNull final SamReader reader, @NotNull final QueryInterval[] intervals) {
        return reader.queryOverlapping(intervals);
    }

    @NotNull
    private static CloseableIterator<SAMRecord> getIterator(@NotNull final SamReader reader, @NotNull final QueryInterval[] intervals,
            @NotNull final long[] filePointers) {
        if (reader instanceof SamReader.PrimitiveSamReaderToSamReaderAdapter) {
            final SamReader.PrimitiveSamReaderToSamReaderAdapter adapter = (SamReader.PrimitiveSamReaderToSamReaderAdapter) reader;
            if (adapter.underlyingReader() instanceof BAMFileReader) {
                final BAMFileReader bamReader = (BAMFileReader) adapter.underlyingReader();
                return bamReader.createIndexIterator(intervals, false, filePointers);
            }
        }
        return reader.queryOverlapping(intervals);
    }

    private static void writeToSlice(@NotNull final String path, @NotNull final SAMFileHeader header,
            @NotNull final CloseableIterator<SAMRecord> iterator) {
        final File outputBAM = new File(path);
        final SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, outputBAM);
        String contig = "";
        while (iterator.hasNext()) {
            final SAMRecord record = iterator.next();
            if (record.getContig() != null && !contig.equals(record.getContig())) {
                contig = record.getContig();
                LOGGER.info("Reading contig: {}", contig);
            }
            writer.addAlignment(record);
        }
        iterator.close();
        writer.close();
    }

    private static int readMaxBufferSize(@NotNull final CommandLine cmd) {
        final String optionValue = cmd.getOptionValue(MAX_CHUNKS_IN_MEMORY);
        try {
            final int bufferSize = Integer.parseInt(optionValue);
            if (bufferSize <= 0) {
                throw new IllegalArgumentException("Buffer size cannot be <= 0.");
            }
            return bufferSize;
        } catch (final NumberFormatException e) {
            throw new IllegalArgumentException("Could not parse buffer size");
        }
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        final OptionGroup inputModeOptionGroup = new OptionGroup();
        inputModeOptionGroup.addOption(Option.builder(INPUT_MODE_S3).required().desc("read input BAM from s3").build());
        inputModeOptionGroup.addOption(Option.builder(INPUT_MODE_FILE).required().desc("read input BAM from file").build());
        inputModeOptionGroup.addOption(Option.builder(INPUT_MODE_URL).required().desc("read input BAM from url").build());
        options.addOptionGroup(inputModeOptionGroup);
        return options;
    }

    @NotNull
    private static Options createURLOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(INPUT_MODE_URL).required().desc("read input BAM from url").build());
        options.addOption(Option.builder(INPUT).required().hasArg().desc("url of BAM file (required)").build());
        options.addOption(Option.builder(INDEX).required().hasArg().desc("url of BAM index file(required)").build());
        return addHttpSlicerOptions(options);
    }

    @NotNull
    private static Options createS3Options() {
        final Options options = new Options();
        options.addOption(Option.builder(INPUT_MODE_S3).required().desc("read input BAM from s3").build());
        options.addOption(Option.builder(BUCKET).required().hasArg().desc("s3 bucket for BAM and index files (required)").build());
        options.addOption(Option.builder(INPUT).required().hasArg().desc("s3 BAM file location (required)").build());
        options.addOption(Option.builder(INDEX).required().hasArg().desc("s3 BAM index location (required)").build());
        return addHttpSlicerOptions(options);
    }

    @NotNull
    private static Options addHttpSlicerOptions(@NotNull final Options options) {
        options.addOption(Option.builder(OUTPUT).required().hasArg().desc("the output BAM (required)").build());
        options.addOption(Option.builder(BED).required().hasArg().desc("BED to slice BAM with (required)").build());
        options.addOption(Option.builder(MAX_CHUNKS_IN_MEMORY).required().hasArg().desc("Max number of chunks to keep in memory").build());
        options.addOption(Option.builder(MAX_CONCURRENT_REQUESTS).required().hasArg().desc("Max concurrent http requests").build());
        return options;
    }

    @NotNull
    private static Options createVcfOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(INPUT_MODE_FILE).required().desc("read input BAM from the filesystem").build());
        options.addOption(Option.builder(INPUT).required().hasArg().desc("the input BAM to slice (required)").build());
        options.addOption(Option.builder(OUTPUT).required().hasArg().desc("the output BAM (required)").build());
        options.addOption(Option.builder(PROXIMITY).hasArg().desc("distance to slice around breakpoint (optional, default=500)").build());
        options.addOption(Option.builder(VCF).required().hasArg().desc("VCF to slice BAM with (required)").build());
        return options;
    }

    @Nullable
    private static CommandLine createCommandLine(@NotNull final String... args) throws ParseException {
        final Options options = createOptions();
        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args, true);
        if (cmd.hasOption(INPUT_MODE_S3)) {
            final Options s3Options = createS3Options();
            try {
                return parser.parse(s3Options, args);

            } catch (ParseException e) {
                LOGGER.error(e.getMessage());
                printHelpAndExit("Slice an s3 BAM file based on BED", s3Options);
            }
        } else if (cmd.hasOption(INPUT_MODE_FILE)) {
            final Options vcfOptions = createVcfOptions();
            try {
                return parser.parse(vcfOptions, args);
            } catch (final ParseException e) {
                LOGGER.error(e.getMessage());
                printHelpAndExit("Slice a local BAM file based on VCF", vcfOptions);
            }
        } else if (cmd.hasOption(INPUT_MODE_URL)) {
            final Options urlOptions = createURLOptions();
            try {
                return parser.parse(urlOptions, args);

            } catch (ParseException e) {
                LOGGER.error(e.getMessage());
                printHelpAndExit("Slice a BAM file based on URLs", urlOptions);
            }
        } else {
            printHelpAndExit("Slice a BAM", options);
        }
        return null;
    }

    private static void printHelpAndExit(@NotNull final String header, @NotNull final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("bam-slicer", header, options, "", true);
        System.exit(1);
    }
}
