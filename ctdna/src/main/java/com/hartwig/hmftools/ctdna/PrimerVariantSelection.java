package com.hartwig.hmftools.ctdna;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.ctdna.PvConfig.PV_LOGGER;
import static com.hartwig.hmftools.ctdna.PvConfig.createCmdLineOptions;
import static com.hartwig.hmftools.ctdna.VariantUtils.calcGcPercent;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.ChromosomeLengthFactory;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.version.VersionInfo;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class PrimerVariantSelection
{
    private final PvConfig mConfig;

    public PrimerVariantSelection(final CommandLine cmd)
    {
        mConfig = new PvConfig(cmd);
    }

    public void run()
    {
        if(!mConfig.isValid())
            System.exit(1);

        PV_LOGGER.info("sample({}) starting primer variant selection", mConfig.SampleId);

        List<Variant> variants = PointMutation.loadSomatics(mConfig);

        variants.addAll(StructuralVariant.loadStructuralVariants(mConfig));

        RefGenomeInterface refGenome = loadRefGenome(mConfig.RefGenomeFile);

        if(refGenome == null)
        {
            PV_LOGGER.error("failed to load ref genome");
            System.exit(1);
        }

        variants.forEach(x -> x.generateSequences(refGenome, mConfig));

        List<Variant> selectedVariants = VariantSelection.selectVariants(variants, mConfig);

        writeSelectedVariants(selectedVariants);

        PV_LOGGER.info("Primer variation selection complete");
    }

    private void writeSelectedVariants(final List<Variant> selectedVariants)
    {
        try
        {
            String filename = mConfig.OutputDir + mConfig.SampleId + ".primer_variants.csv";
            BufferedWriter writer = createBufferedWriter(filename, false);

            writer.write("Category,Variant,CopyNumber,Vaf,TumorFrags,PhasedVariants,Gene");
            writer.write(",Type,Sequence,GcPercent");
            writer.newLine();

            for(Variant variant : selectedVariants)
            {
                String variantInfo = format("%s,%s,%.2f,%.2f,%d,%s,%s",
                        variant.categoryType(), variant.description(), variant.copyNumber(), variant.vaf(),
                        variant.tumorFragments(), variant.hasPhaseVariants(), variant.gene());

                writer.write(format("%s,%s,%s,%.2f", variantInfo, "ALT", variant.sequence(), calcGcPercent(variant.sequence())));
                writer.newLine();

                for(String refSequence : variant.refSequences())
                {
                    writer.write(format("%s,%s,%s,%.2f", variantInfo, "REF", refSequence, calcGcPercent(refSequence)));
                    writer.newLine();
                }
            }

            writer.close();
        }
        catch(IOException e)
        {
            PV_LOGGER.error(" failed to write variants: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        final VersionInfo version = new VersionInfo("ctdna.version");
        PV_LOGGER.info("PrimerVariantSelection version: {}", version.version());

        final Options options = createCmdLineOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            PrimerVariantSelection application = new PrimerVariantSelection(cmd);
            application.run();
        }
        catch(ParseException e)
        {
            PV_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("PrimerVariantSelection", options);
            System.exit(1);
        }
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
