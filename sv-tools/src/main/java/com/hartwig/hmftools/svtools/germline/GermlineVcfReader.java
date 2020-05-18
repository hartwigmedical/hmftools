package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.svtools.germline.AssemblyData.annotateAssembledLinks;
import static com.hartwig.hmftools.svtools.germline.AssemblyData.populateAssemblyLinks;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.AS;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.BEID;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.BEIDL;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.CAS;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.HOMSEQ;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.QUAL;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.RAS;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.REF;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.RP;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.RPQ;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.SR;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.SRQ;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.VF;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.getDoubleValue;
import static com.hartwig.hmftools.svtools.germline.GermlineFilters.getIntValue;
import static com.hartwig.hmftools.svtools.germline.GermlineSV.FIELD_COUNT;
import static com.hartwig.hmftools.svtools.germline.GermlineSV.stripBam;
import static com.hartwig.hmftools.svtools.germline.GermlineVcfConfig.LOG_DEBUG;
import static com.hartwig.hmftools.svtools.germline.GermlineVcfConfig.loadVcfFiles;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;
import com.hartwig.hmftools.common.fusion.GeneAnnotation;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.CommonInfo;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class GermlineVcfReader
{
    private BufferedWriter mCsvWriter;

    private final GermlineFilters mFilter;
    private StructuralVariantFactory mSvFactory;
    private final GermlineVcfConfig mConfig;
    private final GeneImpact mGeneImpact;

    private final List<String> mVcfFiles;

    private final List<GermlineSV> mSampleGermlineSVs;
    private final List<AssemblyData> mSvAssemblyData;

    private static final Logger LOGGER = LogManager.getLogger(GermlineVcfReader.class);

    public GermlineVcfReader(final CommandLine cmd)
    {
        mConfig = new GermlineVcfConfig(cmd);
        mFilter = new GermlineFilters(mConfig);
        mGeneImpact = new GeneImpact(mConfig, cmd);
        mSvFactory = null;

        mVcfFiles = Lists.newArrayList();

        mSvAssemblyData = Lists.newArrayList();
        mSampleGermlineSVs = Lists.newArrayList();
    }

    public void run()
    {
        if(!mConfig.VcfFile.isEmpty())
        {
            mVcfFiles.add(mConfig.VcfFile);
        }
        else if(!mConfig.VcfsFile.isEmpty())
        {
            mVcfFiles.addAll(loadVcfFiles(mConfig.VcfsFile));
        }
        else if(!mConfig.ProcessedFile.isEmpty())
        {
            reprocessVariantsFromFile(mConfig.ProcessedFile);
        }
        else if(!mConfig.BatchRunRootDir.isEmpty())
        {
            findVcfFiles();
        }
        else
        {
            LOGGER.error("missing VCF or batch-run directory");
            return;
        }

        for(final String vcfFile : mVcfFiles)
        {
            processVcf(vcfFile);

            mSvAssemblyData.clear();
        }

        closeBufferedWriter(mCsvWriter);
        mGeneImpact.close();
    }

    private void findVcfFiles()
    {
        // current prod examples
        // structuralVariants/gridss/CPCT02030278R_CPCT02030278T/CPCT02030278R_CPCT02030278T.gridss.vcf.gz
        // structural_caller/WIDE01010356T.gridss.unfiltered.vcf.gz
        final List<String> vcfFiles = Lists.newArrayList();

        try
        {
            final Stream<Path> stream = Files.walk(Paths.get(mConfig.BatchRunRootDir), 5, FileVisitOption.FOLLOW_LINKS);

            mVcfFiles.addAll(stream.filter(x -> !x.toFile().isDirectory())
                    .map(x -> x.toFile().toString())
                    .filter(x -> matchesGridssVcf(x))
                    .collect(Collectors.toList()));
        }
        catch (Exception e)
        {
            LOGGER.error("failed find directories for batchDir({}) run: {}", mConfig.BatchRunRootDir, e.toString());
        }

        LOGGER.info("found {} VCF files", mVcfFiles.size());
    }

    private static boolean matchesGridssVcf(final String filename)
    {
        return filename.endsWith(".gridss.vcf") || filename.endsWith(".gridss.unfiltered.vcf")
                || filename.endsWith(".gridss.vcf.gz") || filename.endsWith(".gridss.unfiltered.vcf.gz");
    }

    private void processVcf(final String vcfFile)
    {
        try
        {
            LOGGER.info("processing germline VCF({})", vcfFile);

            mSampleGermlineSVs.clear();
            mSvFactory = new StructuralVariantFactory(new AlwaysPassFilter());

            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    vcfFile, new VCFCodec(), false);

            reader.iterator().forEach(x -> processVariant(x));

            if(mSampleGermlineSVs.isEmpty())
                return;

            if (mConfig.LinkByAssembly)
                annotateAssembledLinks(mSvAssemblyData);

            if(mConfig.CheckDisruptions)
            {
                final String sampleId = mSampleGermlineSVs.get(0).SampleId;
                mGeneImpact.findDisruptiveVariants(sampleId, mSampleGermlineSVs);
            }

            writeSVs();
        }
        catch(IOException e)
        {
            LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }
    }

    private void processVariant(final VariantContext variant)
    {
        LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

        // early exits
        if(getDoubleValue(variant.getGenotype(0), QUAL) == 0) // no support in the normal
            return;

        if(mConfig.RequirePass && !mConfig.LogFiltered)
        {
            if (!variant.getFilters().isEmpty() && !variant.getFilters().contains(PASS))
                return;
        }

        int currentSvCount = mSvFactory.results().size();
        mSvFactory.addVariantContext(variant);

        // wait for both breakends to be added
        if(currentSvCount == mSvFactory.results().size())
            return;

        // check additional filters
        final StructuralVariant sv = getLastSv();

        if(sv == null)
            return;

        // optionally filter out all but specified chromosomes
        if(!mConfig.RestrictedChromosomes.isEmpty()
        && !mConfig.RestrictedChromosomes.contains(sv.chromosome(true))
        && (sv.type() == SGL || !mConfig.RestrictedChromosomes.contains(sv.chromosome(false))))
        {
            purgeLastSv();
            return;
        }

        List<List<GeneAnnotation>> breakendPairGenes = Lists.newArrayListWithExpectedSize(2);
        breakendPairGenes.add(Lists.newArrayList());
        breakendPairGenes.add(Lists.newArrayList());

        List<GeneAnnotation> overlapGenes = Lists.newArrayList();
        mGeneImpact.populateGeneAnnotations(sv, currentSvCount, breakendPairGenes, overlapGenes);

        if(mConfig.RequireGene && breakendPairGenes.get(SE_START).isEmpty() && breakendPairGenes.get(SE_END).isEmpty() && overlapGenes.isEmpty())
        {
            purgeLastSv();
            return;
        }

        final String filterStr = applyFilters(sv, variant);

        if(!filterStr.equals(PASS) && !mConfig.LogFiltered)
        {
            purgeLastSv();
            return;
        }

        final Genotype normalGenotype = variant.getGenotype(0);
        final String sampleName = stripBam(normalGenotype.getSampleName());

        boolean sglBreakend = sv.type() == SGL;
        final String homology = variant.getAttributeAsString(HOMSEQ, "");

        GermlineSV germlineSV = new GermlineSV(
                sampleName, sv.id(), filterStr, sv.filter(), getDoubleValue(normalGenotype, QUAL), sv.type(),
                sv.chromosome(true), !sglBreakend ? sv.chromosome(false) : "0",
                sv.position(true), !sglBreakend ? sv.position(false) : -1,
                sv.orientation(true), !sglBreakend ? sv.orientation(false) : 0,
                getIntValue(normalGenotype, REF), getIntValue(normalGenotype, RP), getDoubleValue(normalGenotype, RPQ),
                getIntValue(normalGenotype, SR), getDoubleValue(normalGenotype, SRQ), getIntValue(normalGenotype, VF),
                sv.insertSequence(), homology, "", "", "");

        germlineSV.getBreakendGenes().addAll(breakendPairGenes);
        germlineSV.getOverlappedGenes().addAll(overlapGenes);

        mSampleGermlineSVs.add(germlineSV);

        if(mConfig.LinkByAssembly)
            cacheAssemblyData(sv);
    }

    private void cacheAssemblyData(final StructuralVariant sv)
    {
        final CommonInfo variantCI = sv.startContext().getCommonInfo();

        if(variantCI.getAttributeAsInt(AS, 0) + variantCI.getAttributeAsInt(CAS, 0)
        + variantCI.getAttributeAsInt(RAS, 0) < 2)
        {
            return;
        }

        // cache assembly info
        final String startBeId = sv.startContext().getAttributeAsString(BEID,"");
        final String startBeIdl = sv.startContext().getAttributeAsString(BEIDL,"");

        final String endBeId = sv.endContext() != null ? sv.endContext().getAttributeAsString(BEID,"") : "";
        final String endBeIdl = sv.endContext() != null ? sv.endContext().getAttributeAsString(BEIDL,"") : "";

        if(startBeId.isEmpty() && endBeId.isEmpty())
            return;

        final String[] beIdStr = {startBeId, endBeId};
        final String[] beIdlStr = {startBeIdl, endBeIdl};

        AssemblyData asmData = new AssemblyData(sv.id(), beIdStr, beIdlStr);
        mSvAssemblyData.add(asmData);
    }

    private String applyFilters(final StructuralVariant sv, final VariantContext variant)
    {
        if(!variant.getFilters().isEmpty() && mConfig.RequirePass)
            return "GRIDSS_FILTERED";

        if(GermlineFilters.isImprecise(variant))
            return "IMPRECISE";

        // disabled since always zero - see comments from Daniel
        // if(GermlineFilters.noASRP(variant))
        //    return "NO_ASRP";

        if(GermlineFilters.invalidAF(variant))
            return "LOW_AF";

        if(mFilter.belowQualScoreThreshold(variant))
            return "LOW_QS";

        if(GermlineFilters.zeroDiscordantReadSupport(sv, variant))
            return "NO_DISC_RS";

        if(GermlineFilters.zeroSplitReadSupport(sv, variant))
            return "NO_SPLIT_RS";

        if(GermlineFilters.hasStrandBias(sv, variant))
            return "STRAND_BIAS";

        if(GermlineFilters.longPolyCorG(sv, variant))
            return "LONG_POLY_GC";

        return PASS;
    }

    private void writeSVs()
    {
        for(final GermlineSV germlineSV : mSampleGermlineSVs)
        {
            populateAssemblyLinks(mSvAssemblyData, germlineSV);
            writeCsv(germlineSV);
        }
    }

    private void writeCsv(GermlineSV germlineSV)
    {
        try
        {
            if(mCsvWriter == null)
            {
                String outputFileName = !mConfig.VcfFile.isEmpty() ?
                        GermlineSV.generateFilename(mConfig.OutputDir, mConfig.SampleId)
                        : mConfig.OutputDir + "LNX_GERMLINE_SVS.csv";

                mCsvWriter = createBufferedWriter(outputFileName, false);
                mCsvWriter.write(GermlineSV.header());

                mCsvWriter.newLine();
            }

            mCsvWriter.write(GermlineSV.toString(germlineSV));

            mCsvWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing CSV output file: {}", e.toString());
        }
    }

    private void reprocessVariantsFromFile(final String processedFile)
    {
        // may be temporary - read in germline SVs and check for gene disruptions,
        // then write results back out to file
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(processedFile));

            String currentSampleId = "";
            List<GermlineSV> germlineSVs = Lists.newArrayList();

            String line = fileReader.readLine();
            line = fileReader.readLine(); // skip header

            while (line != null)
            {
                final String[] items = line.split(",", -1);

                if(items.length < FIELD_COUNT)
                {
                    LOGGER.error("invalid items: {}", line);
                    return;
                }

                final String sampleId = stripBam(items[0]);

                if(!sampleId.equals(currentSampleId))
                {
                    if(!germlineSVs.isEmpty())
                    {
                        reprocessVariants(currentSampleId, germlineSVs);
                        germlineSVs = Lists.newArrayList();
                    }

                    currentSampleId = sampleId;
                }

                try
                {
                    germlineSVs.add(GermlineSV.fromString(line, true));
                }
                catch(Exception e)
                {
                    LOGGER.error("sample({}) failed to parse line: {}", sampleId, line);
                }

                line = fileReader.readLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error reading file({}): {}", processedFile, e.toString());
        }
    }

    private void reprocessVariants(final String sampleId, final List<GermlineSV> germlineSVs)
    {
        LOGGER.info("sample({}) reprocessing {} germline SVs for disruptions", sampleId, germlineSVs.size());

        int svIndex = 0;
        for(final GermlineSV germlineSV : germlineSVs)
        {
            germlineSV.createSV();
            mGeneImpact.populateGeneAnnotations(germlineSV.sv(), svIndex++, germlineSV.getBreakendGenes(), germlineSV.getOverlappedGenes());
        }

        mGeneImpact.findDisruptiveVariants(sampleId, germlineSVs);

        int disruptiveSvCount = (int)germlineSVs.stream().filter(x -> !x.getDisruptions().isEmpty()).count();

        if(disruptiveSvCount > 0)
        {
            LOGGER.debug("sample({}) has {} disruptive SVs from total({})",
                    germlineSVs.get(0).SampleId, disruptiveSvCount, germlineSVs.size());
        }

        germlineSVs.forEach(x -> writeCsv(x));
    }

    private final StructuralVariant getLastSv()
    {
        return !mSvFactory.results().isEmpty() ? mSvFactory.results().get(mSvFactory.results().size() - 1) : null;
    }

    private void purgeLastSv()
    {
        mSvFactory.results().remove(mSvFactory.results().size() - 1);
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        GermlineVcfConfig.addCommandLineOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        GermlineVcfReader gridssVcfFilters = new GermlineVcfReader(cmd);

        gridssVcfFilters.run();

        LOGGER.info("VCF processing complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
