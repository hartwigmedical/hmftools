package com.hartwig.hmftools.bamslicer;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class BamSlicerApplication {

    private static final String INPUT = "input";
    private static final String OUTPUT = "output";
    private static final String PROXIMITY = "proximity";
    private static final String VCF = "vcf";

    private static Options createOptions() {
        final Options options = new Options();
        options.addOption(Option.builder(INPUT).required().hasArg().desc("the input BAM to slice (required)").build());
        options.addOption(Option.builder(OUTPUT).required().hasArg().desc("the output BAM (required)").build());
        options.addOption(Option.builder(VCF).required().hasArg().desc("VCF to slice BAM with (required)").build());
        options.addOption(Option.builder(PROXIMITY).hasArg().desc("distance to slice around breakpoint (optional, default=500)").build());
        return options;
    }

    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

    private static void printHelpAndExit(final Options options) {
        final HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("Bam-Slicer", "Slice a BAM from a VCF", options, "", true);
        System.exit(1);
    }

    public static void main(final String... args) throws ParseException, IOException {

        final Options options = createOptions();
        try {
            final CommandLine cmd = createCommandLine(options, args);

            // grab arguments
            final String inputPath = cmd.getOptionValue(INPUT);
            final String outputPath = cmd.getOptionValue(OUTPUT);
            final String vcfPath = cmd.getOptionValue(VCF);
            final int proximity = Integer.parseInt(cmd.getOptionValue(PROXIMITY, "500"));

            if (inputPath == null || outputPath == null || vcfPath == null) {
                printHelpAndExit(options);
                return;
            }

            // load the files
            final File inputBAM = new File(inputPath);
            final SamReader reader = SamReaderFactory.makeDefault().open(inputBAM);

            final File vcfFile = new File(vcfPath);
            final VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);

            final List<QueryInterval> queryIntervals = Lists.newArrayList();

            for (VariantContext variant : vcfReader) {

                queryIntervals.add(new QueryInterval(reader.getFileHeader().getSequenceIndex(variant.getContig()),
                        Math.max(0, variant.getStart() - proximity), variant.getStart() + proximity));

                if (variant.getStructuralVariantType() == StructuralVariantType.BND) {

                    final String call = variant.getAlternateAllele(0).getDisplayString();
                    final String[] leftSplit = call.split("\\]");
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
                        System.err.println(variant.getID() + " : could not parse breakpoint");
                        continue;
                    }

                    queryIntervals.add(new QueryInterval(reader.getFileHeader().getSequenceIndex(contig), Math.max(0, position - proximity),
                            position + proximity));
                } else {
                    queryIntervals.add(new QueryInterval(reader.getFileHeader().getSequenceIndex(variant.getContig()),
                            Math.max(0, variant.getEnd() - proximity), variant.getEnd() + proximity));
                }
            }

            final QueryInterval[] optimizedIntervals =
                    QueryInterval.optimizeIntervals(queryIntervals.toArray(new QueryInterval[queryIntervals.size()]));

            writeToSlice(outputPath, reader, optimizedIntervals);
            reader.close();

        } catch (ParseException e) {
            printHelpAndExit(options);
            System.exit(1);
        }
    }

    private static void writeToSlice(final String path, final SamReader reader, final QueryInterval[] intervals) {
        final File outputBAM = new File(path);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), true, outputBAM);
        final SAMRecordIterator iterator = reader.queryOverlapping(intervals);
        while (iterator.hasNext()) {
            writer.addAlignment(iterator.next());
        }
        iterator.close();
        writer.close();
    }

}
