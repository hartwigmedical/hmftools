package com.hartwig.hmftools.bamslicer;

import static htsjdk.samtools.util.BlockCompressedFilePointerUtil.MAX_BLOCK_ADDRESS;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.List;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.slicing.Slicer;
import com.hartwig.hmftools.common.genome.slicing.SlicerFactory;

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
import htsjdk.samtools.cram.CRAIIndex;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.BlockCompressedFilePointerUtil;
import htsjdk.samtools.util.BlockCompressedStreamConstants;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.HttpUtils;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import okhttp3.OkHttpClient;

public class BamSlicerApplication {

    private static final Logger LOGGER = LogManager.getLogger(BamSlicerApplication.class);

    private static final int S3_EXPIRATION_HOURS = 2;

    private static final String INPUT_MODE_S3 = "s3";
    private static final String INPUT_MODE_URL = "url";
    private static final String INPUT_MODE_FILE = "file";

    private static final String S3_ENDPOINT_URL = "s3_endpoint";
    private static final String S3_PROFILE = "s3_profile";
    private static final String S3_BUCKET = "bucket";
    private static final String INPUT = "input";
    private static final String INDEX = "index";
    private static final String OUTPUT = "output";
    private static final String PROXIMITY = "proximity";
    private static final String VCF = "vcf";
    private static final String BED = "bed";
    private static final String UNMAPPED = "unmapped";
    private static final String MAX_CHUNKS_IN_MEMORY = "max_chunks";
    private static final String MAX_CHUNKS_IN_MEMORY_DEFAULT = "2000";
    private static final String MAX_CONCURRENT_REQUESTS = "max_concurrent_requests";
    private static final String MAX_CONCURRENT_REQUESTS_DEFAULT = "50";

    private static final String REF_GENOME_FASTA_FILE = "ref_genome_fasta_file";

    private static final Chunk HEADER_CHUNK = new Chunk(0, (long) BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE << 16);

    public static void main(final String... args) throws ParseException, IOException {
        CommandLine cmd = createCommandLine(args);

        // Disable default samtools buffering
        System.setProperty("samjdk.buffer_size", "0");
        if (cmd.hasOption(INPUT_MODE_FILE)) {
            sliceFromVCF(cmd);
        } else if (cmd.hasOption(INPUT_MODE_S3)) {
            Pair<URL, URL> urls = generateURLs(cmd);
            sliceFromURLs(urls.getValue(), urls.getKey(), cmd);
        } else if (cmd.hasOption(INPUT_MODE_URL)) {
            URL bamURL = new URL(cmd.getOptionValue(INPUT));
            URL indexURL = new URL(cmd.getOptionValue(INDEX));
            sliceFromURLs(indexURL, bamURL, cmd);
        }

        LOGGER.info("Done.");
    }

    private static void sliceFromVCF(@NotNull CommandLine cmd) throws IOException {
        String inputPath = cmd.getOptionValue(INPUT);
        String vcfPath = cmd.getOptionValue(VCF);
        int proximity = Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500"));

        SamReaderFactory readerFactory = createFromCommandLine(cmd);
        SamReader reader = readerFactory.open(new File(inputPath));

        QueryInterval[] intervals = getIntervalsFromVCF(vcfPath, reader.getFileHeader(), proximity);
        CloseableIterator<SAMRecord> iterator = reader.queryOverlapping(intervals);
        SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(reader.getFileHeader(), true, new File(cmd.getOptionValue(OUTPUT)));

        writeToSlice(writer, iterator);

        writer.close();
        reader.close();
    }

    @NotNull
    private static QueryInterval[] getIntervalsFromVCF(@NotNull String vcfPath, @NotNull SAMFileHeader header, int proximity) {
        File vcfFile = new File(vcfPath);
        VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);
        List<QueryInterval> queryIntervals = Lists.newArrayList();

        for (VariantContext variant : vcfReader) {
            queryIntervals.add(new QueryInterval(header.getSequenceIndex(variant.getContig()),
                    Math.max(0, variant.getStart() - proximity),
                    variant.getStart() + proximity));

            if (variant.getStructuralVariantType() == StructuralVariantType.BND) {
                String call = variant.getAlternateAllele(0).getDisplayString();
                String[] leftSplit = call.split("]");
                String[] rightSplit = call.split("\\[");

                String contig;
                int position;
                if (leftSplit.length >= 2) {
                    final String[] location = leftSplit[1].split(":");
                    contig = location[0];
                    position = Integer.parseInt(location[1]);
                } else if (rightSplit.length >= 2) {
                    final String[] location = rightSplit[1].split(":");
                    contig = location[0];
                    position = Integer.parseInt(location[1]);
                } else {
                    LOGGER.error("{} : could not parse breakpoint", variant.getID());
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
    private static Pair<URL, URL> generateURLs(@NotNull CommandLine cmd) {
        try {
            String s3Endpoint = cmd.getOptionValue(S3_ENDPOINT_URL);
            String s3Profile = cmd.getOptionValue(S3_PROFILE);

            LOGGER.info("Attempting to generate S3 URLs for endpoint: {} using profile: {}", s3Endpoint, s3Profile);
            S3UrlGenerator urlGenerator = ImmutableS3UrlGenerator.of(s3Endpoint, s3Profile);
            URL bamUrl = urlGenerator.generateUrl(cmd.getOptionValue(S3_BUCKET), cmd.getOptionValue(INPUT), S3_EXPIRATION_HOURS);
            URL indexUrl = urlGenerator.generateUrl(cmd.getOptionValue(S3_BUCKET), cmd.getOptionValue(INDEX), S3_EXPIRATION_HOURS);
            return Pair.of(bamUrl, indexUrl);
        } catch (Exception e) {
            LOGGER.error("Could not create S3 URLs. Error: {}", e.toString());
            System.exit(1);
        }
        return null;
    }

    private static void sliceFromURLs(@NotNull URL indexUrl, @NotNull URL bamUrl, @NotNull CommandLine cmd) throws IOException {
        File indexFile = downloadIndex(indexUrl);
        indexFile.deleteOnExit();

        SamReader reader = createFromCommandLine(cmd).open(SamInputResource.of(bamUrl).index(indexFile));

        BAMIndex bamIndex;
        if (indexFile.getPath().contains(".crai")) {
            SeekableStream craiIndex = CRAIIndex.openCraiFileAsBaiStream(indexFile, reader.getFileHeader().getSequenceDictionary());
            bamIndex = new DiskBasedBAMFileIndex(craiIndex, reader.getFileHeader().getSequenceDictionary());
        } else {
            bamIndex = new DiskBasedBAMFileIndex(indexFile, reader.getFileHeader().getSequenceDictionary(), false);
        }

        Optional<Pair<QueryInterval[], BAMFileSpan>> queryIntervalsAndSpan = queryIntervalsAndSpan(reader, bamIndex, cmd);
        Optional<Chunk> unmappedChunk = getUnmappedChunk(bamIndex, HttpUtils.getHeaderField(bamUrl, "Content-Length"), cmd);
        List<Chunk> sliceChunks = sliceChunks(queryIntervalsAndSpan, unmappedChunk);
        SamReader cachingReader = createCachingReader(indexFile, bamUrl, cmd, sliceChunks);

        SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true)
                .makeBAMWriter(reader.getFileHeader(), true, new File(cmd.getOptionValue(OUTPUT)));

        queryIntervalsAndSpan.ifPresent(pair -> {
            LOGGER.info("Slicing bam on bed regions...");
            CloseableIterator<SAMRecord> bedIterator = getIterator(cachingReader, pair.getKey(), pair.getValue().toCoordinateArray());
            writeToSlice(writer, bedIterator);
            LOGGER.info("Done writing bed slices.");
        });

        unmappedChunk.ifPresent(chunk -> {
            LOGGER.info("Slicing unmapped reads...");
            CloseableIterator<SAMRecord> unmappedIterator = cachingReader.queryUnmapped();
            writeToSlice(writer, unmappedIterator);
            LOGGER.info("Done writing unmapped reads.");
        });

        reader.close();
        writer.close();
        cachingReader.close();
    }

    @NotNull
    private static Optional<Pair<QueryInterval[], BAMFileSpan>> queryIntervalsAndSpan(@NotNull SamReader reader, @NotNull BAMIndex bamIndex,
            @NotNull CommandLine cmd) throws IOException {
        if (cmd.hasOption(BED)) {
            String bedPath = cmd.getOptionValue(BED);
            LOGGER.info("Reading query intervals from BED file: {}", bedPath);
            QueryInterval[] intervals = getIntervalsFromBED(bedPath, reader.getFileHeader());
            BAMFileSpan span = BAMFileReader.getFileSpan(intervals, bamIndex);
            return Optional.of(Pair.of(intervals, span));
        }
        return Optional.empty();
    }

    @NotNull
    private static SamReader createCachingReader(@NotNull File indexFile, @NotNull URL bamUrl, @NotNull CommandLine cmd,
            @NotNull List<Chunk> sliceChunks) throws IOException {
        OkHttpClient httpClient =
                SlicerHttpClient.create(Integer.parseInt(cmd.getOptionValue(MAX_CONCURRENT_REQUESTS, MAX_CONCURRENT_REQUESTS_DEFAULT)));
        int maxBufferSize = readMaxBufferSize(cmd);

        SamInputResource bamResource =
                SamInputResource.of(new CachingSeekableHTTPStream(httpClient, bamUrl, sliceChunks, maxBufferSize)).index(indexFile);
        SamReaderFactory readerFactory = createFromCommandLine(cmd);

        return readerFactory.open(bamResource);
    }

    @NotNull
    private static List<Chunk> sliceChunks(@NotNull Optional<Pair<QueryInterval[], BAMFileSpan>> queryIntervalsAndSpan,
            @NotNull Optional<Chunk> unmappedChunk) {
        List<Chunk> chunks = Lists.newArrayList();
        chunks.add(HEADER_CHUNK);
        queryIntervalsAndSpan.ifPresent(pair -> {
            chunks.addAll(expandChunks(pair.getValue().getChunks()));
            LOGGER.info("Generated {} query intervals which map to {} bam chunks", pair.getKey().length, chunks.size());
        });
        unmappedChunk.ifPresent(chunks::add);
        return Chunk.optimizeChunkList(chunks, 0);
    }

    @NotNull
    private static List<Chunk> expandChunks(@NotNull List<Chunk> chunks) {
        List<Chunk> result = Lists.newArrayList();
        for (Chunk chunk : chunks) {
            long chunkEndBlockAddress = BlockCompressedFilePointerUtil.getBlockAddress(chunk.getChunkEnd());
            long extendedEndBlockAddress = chunkEndBlockAddress + BlockCompressedStreamConstants.MAX_COMPRESSED_BLOCK_SIZE;
            long newChunkEnd = Math.min(extendedEndBlockAddress, MAX_BLOCK_ADDRESS);
            long chunkEndVirtualPointer = newChunkEnd << 16;
            result.add(new Chunk(chunk.getChunkStart(), chunkEndVirtualPointer));
        }
        return result;
    }

    @NotNull
    private static File downloadIndex(@NotNull URL indexUrl) throws IOException {
        LOGGER.info("Downloading index from {}", indexUrl);

        ReadableByteChannel indexChannel = Channels.newChannel(indexUrl.openStream());

        String extension = ".bai";
        if (indexUrl.getPath().contains(".crai")) {
            LOGGER.debug("Located crai from the index on {}", indexUrl);
            extension = ".crai";
        }

        File index = File.createTempFile("tmp", extension);
        FileOutputStream indexOutputStream = new FileOutputStream(index);
        indexOutputStream.getChannel().transferFrom(indexChannel, 0, Long.MAX_VALUE);
        indexOutputStream.close();
        indexChannel.close();

        LOGGER.info("Downloaded index to {}", index.getPath());

        return index;
    }

    @NotNull
    private static Optional<Chunk> getUnmappedChunk(@NotNull BAMIndex bamIndex, @Nullable String contentLengthString,
            @NotNull CommandLine cmd) {
        if (cmd.hasOption(UNMAPPED)) {
            long startOfLastLinearBin = bamIndex.getStartOfLastLinearBin();
            if (startOfLastLinearBin == -1) {
                LOGGER.warn("Start of last linear bin was -1. No mapped reads found in BAM.");
                return Optional.empty();
            }
            if (contentLengthString != null) {
                try {
                    long contentLength = Long.parseLong(contentLengthString);
                    // We multiply content length with 2^16 = ~64k. Presumably 'content length' is "in terms of number of 64Kb packets".
                    return Optional.of(new Chunk(startOfLastLinearBin, contentLength << 16));
                } catch (NumberFormatException ignored) {
                    LOGGER.error("Invalid content length ({}) for bam URL", contentLengthString);
                    return Optional.empty();
                }
            }
        }
        return Optional.empty();
    }

    @NotNull
    private static QueryInterval[] getIntervalsFromBED(@NotNull String bedPath, @NotNull SAMFileHeader header) throws IOException {
        Slicer bedSlicer = SlicerFactory.fromBedFile(bedPath);
        List<QueryInterval> queryIntervals = Lists.newArrayList();
        for (GenomeRegion region : bedSlicer.regions()) {
            queryIntervals.add(new QueryInterval(header.getSequenceIndex(region.chromosome()), (int) region.start(), (int) region.end()));
        }
        return QueryInterval.optimizeIntervals(queryIntervals.toArray(new QueryInterval[queryIntervals.size()]));
    }

    @NotNull
    private static CloseableIterator<SAMRecord> getIterator(@NotNull SamReader reader, @NotNull QueryInterval[] intervals,
            @NotNull final long[] filePointers) {
        if (reader instanceof SamReader.PrimitiveSamReaderToSamReaderAdapter) {
            SamReader.PrimitiveSamReaderToSamReaderAdapter adapter = (SamReader.PrimitiveSamReaderToSamReaderAdapter) reader;
            if (adapter.underlyingReader() instanceof BAMFileReader) {
                BAMFileReader bamReader = (BAMFileReader) adapter.underlyingReader();
                return bamReader.createIndexIterator(intervals, false, filePointers);
            }
        }
        return reader.queryOverlapping(intervals);
    }

    private static void writeToSlice(@NotNull SAMFileWriter writer, @NotNull CloseableIterator<SAMRecord> iterator) {
        String contig = "";
        while (iterator.hasNext()) {
            SAMRecord record = iterator.next();
            if (record.getContig() != null && !contig.equals(record.getContig())) {
                contig = record.getContig();
                LOGGER.info("Reading contig: {}", contig);
            }
            writer.addAlignment(record);
        }
        iterator.close();
    }

    private static int readMaxBufferSize(@NotNull final CommandLine cmd) {
        String optionValue = cmd.getOptionValue(MAX_CHUNKS_IN_MEMORY, MAX_CHUNKS_IN_MEMORY_DEFAULT);
        try {
            int bufferSize = Integer.parseInt(optionValue);
            if (bufferSize <= 0) {
                throw new IllegalArgumentException("Buffer size cannot be <= 0.");
            }
            return bufferSize;
        } catch (final NumberFormatException e) {
            throw new IllegalArgumentException("Could not parse buffer size");
        }
    }

    @NotNull
    private static SamReaderFactory createFromCommandLine(@NotNull CommandLine cmd) {
        SamReaderFactory readerFactory = SamReaderFactory.makeDefault();

        if (cmd.hasOption(REF_GENOME_FASTA_FILE)) {
            readerFactory.referenceSource(new ReferenceSource(new File(cmd.getOptionValue(REF_GENOME_FASTA_FILE))));
        }

        return readerFactory;
    }

    @NotNull
    private static Options createOptions() {
        Options options = new Options();
        OptionGroup inputModeOptionGroup = new OptionGroup();
        inputModeOptionGroup.addOption(Option.builder(INPUT_MODE_S3).required().desc("read input BAM from s3").build());
        inputModeOptionGroup.addOption(Option.builder(INPUT_MODE_FILE).required().desc("read input BAM from file").build());
        inputModeOptionGroup.addOption(Option.builder(INPUT_MODE_URL).required().desc("read input BAM from url").build());
        options.addOptionGroup(inputModeOptionGroup);
        return options;
    }

    @NotNull
    private static Options createURLOptions() {
        Options options = new Options();
        options.addOption(Option.builder(INPUT_MODE_URL).required().desc("Read input BAM from url").build());
        options.addOption(Option.builder(INPUT).required().hasArg().desc("url of BAM file (required)").build());
        options.addOption(Option.builder(INDEX).required().hasArg().desc("url of BAM index file(required)").build());
        options.addOption(Option.builder(REF_GENOME_FASTA_FILE)
                .hasArg()
                .desc("(Optional) path to the ref genome fasta (for reading CRAMs)")
                .build());
        return addHttpSlicerOptions(options);
    }

    @NotNull
    private static Options createS3Options() {
        Options options = new Options();
        options.addOption(Option.builder(INPUT_MODE_S3).required().desc("Read input BAM from s3").build());
        options.addOption(Option.builder(S3_ENDPOINT_URL).required().hasArg().desc("URL of the s3 end point").build());
        options.addOption(Option.builder(S3_PROFILE).required().hasArg().desc("Name of the profile to access the s3 bucket").build());
        options.addOption(Option.builder(S3_BUCKET).required().hasArg().desc("s3 bucket for BAM and index files (required)").build());
        options.addOption(Option.builder(INPUT).required().hasArg().desc("s3 BAM file location (required)").build());
        options.addOption(Option.builder(INDEX).required().hasArg().desc("s3 BAM index location (required)").build());
        return addHttpSlicerOptions(options);
    }

    @NotNull
    private static Options addHttpSlicerOptions(@NotNull final Options options) {
        options.addOption(Option.builder(OUTPUT).required().hasArg().desc("The output BAM (required)").build());
        options.addOption(Option.builder(BED).hasArg().desc("BED to slice BAM with").build());
        options.addOption(Option.builder(UNMAPPED).desc("Slice unmapped reads").build());
        options.addOption(Option.builder(MAX_CHUNKS_IN_MEMORY)
                .hasArg()
                .desc("Max number of chunks to keep in memory (default: " + MAX_CHUNKS_IN_MEMORY_DEFAULT + ")")
                .build());
        options.addOption(Option.builder(MAX_CONCURRENT_REQUESTS)
                .hasArg()
                .desc("Max concurrent http requests (default: " + MAX_CONCURRENT_REQUESTS_DEFAULT + ")")
                .build());
        return options;
    }

    @NotNull
    private static Options createVcfOptions() {
        Options options = new Options();
        options.addOption(Option.builder(INPUT_MODE_FILE).required().desc("read input BAM from the filesystem").build());
        options.addOption(Option.builder(INPUT).required().hasArg().desc("the input BAM to slice (required)").build());
        options.addOption(Option.builder(OUTPUT).required().hasArg().desc("the output BAM (required)").build());
        options.addOption(Option.builder(PROXIMITY).hasArg().desc("distance to slice around breakpoint (optional, default=500)").build());
        options.addOption(Option.builder(VCF).required().hasArg().desc("VCF to slice BAM with (required)").build());
        options.addOption(Option.builder(REF_GENOME_FASTA_FILE)
                .hasArg()
                .desc("(Optional) path to the ref genome fasta (for reading CRAMs)")
                .build());
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull String... args) throws ParseException {
        Options options = createOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmd = parser.parse(options, args, true);
        if (cmd.hasOption(INPUT_MODE_S3)) {
            Options s3Options = createS3Options();
            try {
                return parser.parse(s3Options, args);
            } catch (ParseException e) {
                LOGGER.error(e.getMessage());
                printHelpAndExit("Slice an s3 BAM file based on BED", s3Options);
            }
        } else if (cmd.hasOption(INPUT_MODE_FILE)) {
            Options vcfOptions = createVcfOptions();
            try {
                return parser.parse(vcfOptions, args);
            } catch (final ParseException e) {
                LOGGER.error(e.getMessage());
                printHelpAndExit("Slice a local BAM file based on VCF", vcfOptions);
            }
        } else if (cmd.hasOption(INPUT_MODE_URL)) {
            Options urlOptions = createURLOptions();
            try {
                return parser.parse(urlOptions, args);
            } catch (ParseException e) {
                LOGGER.error(e.getMessage());
                printHelpAndExit("Slice a BAM file based on URLs", urlOptions);
            }
        } else {
            printHelpAndExit("Slice a BAM", options);
        }

        throw new IllegalStateException("Option creation reached end of options. Shouldn't be possible");
    }

    private static void printHelpAndExit(@NotNull String header, @NotNull Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("bam-slicer", header, options, "", true);
        System.exit(1);
    }
}
