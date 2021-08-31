package com.hartwig.hmftools.svtools.pon;

import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.svtools.pon.PonLocations.chrEnd;
import static com.hartwig.hmftools.svtools.pon.PonLocations.chrStart;
import static com.hartwig.hmftools.svtools.pon.PonLocations.chromosome;
import static com.hartwig.hmftools.svtools.pon.PonLocations.locationValues;
import static com.hartwig.hmftools.svtools.pon.PonLocations.orientEnd;
import static com.hartwig.hmftools.svtools.pon.PonLocations.orientStart;
import static com.hartwig.hmftools.svtools.pon.PonLocations.orientation;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class PonBuilder
{
    private StructuralVariantFactory mSvFactory;

    private final PonConfig mConfig;
    private final Map<String,String> mSampleVcfFiles;
    private final PonStore mPonStore;

    private int mProcessedVariants;

    public static final Logger PON_LOGGER = LogManager.getLogger(PonBuilder.class);

    public PonBuilder(final CommandLine cmd)
    {
        mConfig = new PonConfig(cmd);
        mSvFactory = null;

        mSampleVcfFiles = Maps.newHashMap();
        mPonStore = new PonStore();

        mProcessedVariants = 0;
    }

    public void run()
    {
        PON_LOGGER.info("building PON from {} samples", mConfig.SampleIds.size());

        findVcfFiles();

        if(mSampleVcfFiles.isEmpty())
        {
            PON_LOGGER.error("missing VCF or batch-run directory");
            System.exit(1);
        }

        PON_LOGGER.info("found {} sample VCFs", mSampleVcfFiles.size());

        for(Map.Entry<String,String> entry : mSampleVcfFiles.entrySet())
        {
            processVcf(entry.getKey(), entry.getValue());
        }

        PON_LOGGER.info("writing PON files");

        writePonFiles();

        PON_LOGGER.info("PON creation complete");
    }

    private void findVcfFiles()
    {
        try
        {
            final Stream<Path> stream = Files.walk(Paths.get(mConfig.RootDirectory), 3, FileVisitOption.FOLLOW_LINKS);

            List<String> candidateVcfFiles = stream.filter(x -> !x.toFile().isDirectory())
                    .map(x -> x.toFile().toString())
                    .filter(x -> mConfig.VcfFileIds.stream().anyMatch(y -> x.endsWith(y)))
                    .collect(Collectors.toList());

            for(String sampleId : mConfig.SampleIds)
            {
                String sampleVcfFile = candidateVcfFiles.stream().filter(x -> x.contains(sampleId)).findFirst().orElse(null);

                if(sampleVcfFile != null)
                    mSampleVcfFiles.put(sampleId, sampleVcfFile);
            }

            PON_LOGGER.info("found {} VCF files", mSampleVcfFiles.size());
        }
        catch (Exception e)
        {
            PON_LOGGER.error("failed find directories for batchDir({}) run: {}", mConfig.RootDirectory, e.toString());
        }
    }

    private void clearSampleData()
    {
        mProcessedVariants = 0;
        mSvFactory = new StructuralVariantFactory(new AlwaysPassFilter());
    }

    private void processVcf(final String sampleId, final String vcfFile)
    {
        clearSampleData();

        try
        {
            PON_LOGGER.info("processing sampleId({}) vcfFie({})", sampleId, vcfFile);

            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    vcfFile, new VCFCodec(), false);

            reader.iterator().forEach(x -> processVariant(x));

        }
        catch(IOException e)
        {
            PON_LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }

        PON_LOGGER.info("sample({}) read {} variants", sampleId, mProcessedVariants);

        // log current PON stats
        PON_LOGGER.info("PON stats: SVs(loc={} variants={}) SGLs(loc={} variants={)",
                mPonStore.svLocationCount(), mPonStore.svPonCount(),
                mPonStore.sglLocationCount(), mPonStore.sglPonCount());
    }

    private void processVariant(final VariantContext variant)
    {
        PON_LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

        ++mProcessedVariants;

        if(mProcessedVariants > 0 && (mProcessedVariants % 100000) == 0)
        {
            // PON_LOGGER.debug("sample({}) processed {} variants, VCF-unmatched({})",
            //        mConfig.SampleId, mProcessedVariants, mSvFactory.unmatched().size());
        }

        int currentSvCount = mSvFactory.results().size();
        mSvFactory.addVariantContext(variant);

        // wait for both breakends to be added
        if(currentSvCount == mSvFactory.results().size())
            return;

        final StructuralVariant sv = popLastSv(); // get and clear from storage

        if(sv == null)
            return;

        if(sv.type() == SGL)
        {
            mPonStore.addLocation(sv.chromosome(true), sv.orientation(true), sv.position(true).intValue());
        }
        else
        {
            mPonStore.addLocation(
                    sv.chromosome(true), sv.chromosome(false),
                    sv.orientation(true), sv.orientation(false),
                    sv.position(true).intValue(), sv.position(false).intValue());
        }
    }

    private final StructuralVariant popLastSv()
    {
        if(mSvFactory.results().isEmpty())
            return null;

        StructuralVariant sv = mSvFactory.results().get(0);
        mSvFactory.results().remove(0);

        return sv;
    }

    private void writePonFiles()
    {
        try
        {
            String svPonFile = mConfig.OutputDir + "sv_pon.csv";
            BufferedWriter svWriter = createBufferedWriter(svPonFile, false);
            svWriter.write("ChrStart,OrientStart,PosStart,ChrEnd,OrientEnd,PosEnd,PonCount");
            svWriter.newLine();

            for(Map.Entry<String,PonLocations> entry : mPonStore.getSvLocations().entrySet())
            {
                String locationKey = entry.getKey();
                String[] locValues = locationValues(locationKey);
                PonLocations locations = entry.getValue();

                for(LocationCounter locationStart : locations.Locations)
                {
                    for(LocationCounter locationEnd : locationStart.getNextLocations())
                    {
                        if(locationEnd.getCount() < mConfig.MinPonWriteCount)
                            continue;

                        svWriter.write(String.format("%s,%d,%d,%s,%d,%d,%d",
                                chrStart(locValues), orientStart(locValues), locationStart.Position,
                                chrEnd(locValues), orientEnd(locValues), locationEnd.Position, locationEnd.getCount()));
                        svWriter.newLine();
                    }
                }
            }

            svWriter.close();

            String sglPonFile = mConfig.OutputDir + "sgl_pon.csv";
            BufferedWriter sglWriter = createBufferedWriter(sglPonFile, false);
            sglWriter.write("Chromosome,Orientation,Position,PonCount");
            sglWriter.newLine();

            for(Map.Entry<String,PonLocations> entry : mPonStore.getSglLocations().entrySet())
            {
                String locationKey = entry.getKey();
                String[] locValues = locationValues(locationKey);
                PonLocations locations = entry.getValue();

                for(LocationCounter locationStart : locations.Locations)
                {
                    if(locationStart.getCount() < mConfig.MinPonWriteCount)
                        continue;

                    sglWriter.write(String.format("%s,%d,%d,%d",
                            chromosome(locValues), orientation(locValues), locationStart.Position, locationStart.getCount()));
                    sglWriter.newLine();
                }
            }

            sglWriter.close();
        }
        catch(IOException e)
        {
            PON_LOGGER.error("failed to write PON files: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        PonConfig.addOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PonBuilder germlineVcfReader = new PonBuilder(cmd);
        germlineVcfReader.run();

        PON_LOGGER.info("VCF processing complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
