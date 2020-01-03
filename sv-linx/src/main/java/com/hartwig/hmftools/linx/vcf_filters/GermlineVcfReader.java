package com.hartwig.hmftools.linx.vcf_filters;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.vcf_filters.AssemblyData.annotateAssembledLinks;
import static com.hartwig.hmftools.linx.vcf_filters.AssemblyData.populateAssemblyLinks;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.AS;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.BEID;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.BEIDL;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.CAS;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.HOMSEQ;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.QUAL;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.RAS;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.REF;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.RP;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.RPQ;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.SR;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.SRQ;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.VF;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.getDoubleValue;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.getIntValue;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineSV.FIELD_COUNT;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineVcfConfig.LOG_DEBUG;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineVcfConfig.loadVcfFiles;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.annotation.GeneAnnotation;

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

    private final List<AssemblyData> mSvAssemblyData;
    private final Map<String,String> mSvFilterReasons;

    private static final Logger LOGGER = LogManager.getLogger(GermlineVcfReader.class);

    public GermlineVcfReader(final CommandLine cmd)
    {
        mConfig = new GermlineVcfConfig(cmd);
        mFilter = new GermlineFilters(mConfig);
        mGeneImpact = new GeneImpact(mConfig, cmd);
        mSvFactory = null;

        mVcfFiles = Lists.newArrayList();
        mSvFilterReasons = Maps.newHashMap();
        mSvAssemblyData = Lists.newArrayList();
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
            mSvFilterReasons.clear();
        }

        closeBufferedWriter(mCsvWriter);
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

            mSvFactory = new StructuralVariantFactory(new AlwaysPassFilter());

            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    vcfFile, new VCFCodec(), false);

            reader.iterator().forEach(x -> processVariant(x));

            if (mConfig.LinkByAssembly)
                annotateAssembledLinks(mSvAssemblyData);

            if(mConfig.CheckDisruptions)
                mGeneImpact.findDisruptiveVariants(mSvFactory.results());

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

        if(mConfig.RequireGene && !mGeneImpact.hasGeneAnnotation(sv))
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

        mSvFilterReasons.put(sv.id(), filterStr);

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
        // final Set<StructuralVariant> disruptedSVs = mGeneImpact.getGeneDisruptions().keySet();
        final Map<StructuralVariant,String> disruptionTypes = mGeneImpact.getDisruptionTypes();

        for(final StructuralVariant sv : mSvFactory.results())
        {
            final String filterStr = mSvFilterReasons.get(sv.id());

            final VariantContext variant = sv.startContext();

            final Genotype normalGenotype = variant.getGenotype(0);
            final String sampleName = normalGenotype.getSampleName();
            // final Genotype tumorGenotype = variant.getGenotype(1);

            boolean sglBreakend = sv.type() == SGL;

            final String homology = variant.getAttributeAsString(HOMSEQ, "");

            String[] geneAnnotations = {"",""};
            mGeneImpact.annotationWithGenes(sv, geneAnnotations);

            String[] asmbSvIds = {"",""};
            populateAssemblyLinks(mSvAssemblyData, sv, asmbSvIds);

            String genePanelOverlaps = mGeneImpact.annotateWithGenePanel(sv);

            GermlineSV germlineSV = new GermlineSV(
                    sampleName, sv.id(), filterStr, sv.filter(), getDoubleValue(normalGenotype, QUAL), sv.type(),
                    sv.chromosome(true), !sglBreakend ? sv.chromosome(false) : "0",
                    sv.position(true), !sglBreakend ? sv.position(false) : -1,
                    sv.orientation(true), !sglBreakend ? sv.orientation(false) : 0,
                    getIntValue(normalGenotype, REF), getIntValue(normalGenotype, RP), getDoubleValue(normalGenotype, RPQ),
                    getIntValue(normalGenotype, SR), getDoubleValue(normalGenotype, SRQ), getIntValue(normalGenotype, VF),
                    sv.insertSequence(), homology,
                    geneAnnotations[SE_START], geneAnnotations[SE_END], genePanelOverlaps, asmbSvIds[SE_START], asmbSvIds[SE_END]);

            if(disruptionTypes.containsKey(sv))
                germlineSV.setDisruptionType(disruptionTypes.get(sv));

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

                if(mConfig.CheckDisruptions)
                    mCsvWriter.write(",DisruptionType");

                mCsvWriter.newLine();
            }

            mCsvWriter.write(GermlineSV.toString(germlineSV));

            if(mConfig.CheckDisruptions)
                mCsvWriter.write(String.format(",%s", germlineSV.disruptionType()));

            mCsvWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing CSV output file: {}", e.toString());
        }
    }

    private void reprocessVariantsFromFile(final String processedFile)
    {
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

                final String sampleId = items[0];

                if(!sampleId.equals(currentSampleId))
                {
                    if(!germlineSVs.isEmpty())
                    {
                        reprocessVariants(germlineSVs);
                        germlineSVs = Lists.newArrayList();
                    }

                    currentSampleId = sampleId;
                }

                germlineSVs.add(GermlineSV.fromString(line));
                line = fileReader.readLine();
            }
        }
        catch (final IOException e)
        {
            LOGGER.error("error reading file({}): {}", processedFile, e.toString());
        }
    }

    private void reprocessVariants(final List<GermlineSV> germlineSVs)
    {
        final List<StructuralVariant> svList = Lists.newArrayList();

        germlineSVs.forEach(x -> svList.add(GermlineSV.convert(x)));

        mGeneImpact.findDisruptiveVariants(svList);

        final Map<StructuralVariant,List<GeneAnnotation>> svGeneDisruptions = mGeneImpact.getGeneDisruptions();

        if(!svGeneDisruptions.isEmpty())
        {
            LOGGER.info("sample({}) has {} disruptive SVs from total({})",
                    germlineSVs.get(0).SampleId, svGeneDisruptions.size(), germlineSVs.size());
        }
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
