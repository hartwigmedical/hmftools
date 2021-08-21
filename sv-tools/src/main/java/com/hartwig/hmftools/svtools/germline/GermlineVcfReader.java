package com.hartwig.hmftools.svtools.germline;

import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
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
import static com.hartwig.hmftools.svtools.germline.GermlineUtils.GM_LOGGER;
import static com.hartwig.hmftools.svtools.germline.GermlineUtils.findVcfFiles;
import static com.hartwig.hmftools.svtools.germline.GermlineUtils.stripBam;
import static com.hartwig.hmftools.svtools.germline.GermlineVcfConfig.loadVcfFiles;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
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

    private final List<String> mVcfFiles;
    private final LinkAnalyser mLinkAnalyser;
    private final PonCache mPonCache;

    private final List<GermlineSV> mSampleGermlineSVs;

    public GermlineVcfReader(final CommandLine cmd)
    {
        mConfig = new GermlineVcfConfig(cmd);
        mFilter = new GermlineFilters(mConfig);
        mSvFactory = null;

        mVcfFiles = Lists.newArrayList();
        mLinkAnalyser = new LinkAnalyser();
        mPonCache = new PonCache(cmd);

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
            mVcfFiles.addAll(findVcfFiles(mConfig.BatchRunRootDir));
        }
        else
        {
            GM_LOGGER.error("missing VCF or batch-run directory");
            return;
        }

        for(final String vcfFile : mVcfFiles)
        {
            processVcf(vcfFile);

            mLinkAnalyser.clear();
        }

        closeBufferedWriter(mCsvWriter);
    }

    private void processVcf(final String vcfFile)
    {
        try
        {
            GM_LOGGER.info("processing germline VCF({})", vcfFile);

            mSampleGermlineSVs.clear();
            mSvFactory = new StructuralVariantFactory(new AlwaysPassFilter());

            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    vcfFile, new VCFCodec(), false);

            reader.iterator().forEach(x -> processVariant(x));

            if(mSampleGermlineSVs.isEmpty())
                return;

            if (mConfig.LinkByAssembly)
                mLinkAnalyser.annotateAssembledLinks();

            for(final GermlineSV germlineSV : mSampleGermlineSVs)
            {
                mLinkAnalyser.populateAssemblyLinks(germlineSV);
                writeCsv(germlineSV);
            }
        }
        catch(IOException e)
        {
            GM_LOGGER.error("error reading vcf({}): {}", vcfFile, e.toString());
        }
    }

    private void processVariant(final VariantContext variant)
    {
        GM_LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

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

        final String filterStr = mFilter.applyFilters(sv, variant);

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
        mLinkAnalyser.add(asmData);
    }

    private void writeCsv(GermlineSV germlineSV)
    {
        try
        {
            if(mCsvWriter == null)
            {
                String outputFileName = !mConfig.VcfFile.isEmpty() ?
                        CsvFileWriter.generateFilename(mConfig.OutputDir, mConfig.SampleId)
                        : mConfig.OutputDir + "LNX_GERMLINE_SVS.csv";

                mCsvWriter = createBufferedWriter(outputFileName, false);
                mCsvWriter.write(CsvFileWriter.header());

                mCsvWriter.newLine();
            }

            mCsvWriter.write(CsvFileWriter.toString(germlineSV));

            mCsvWriter.newLine();
        }
        catch (final IOException e)
        {
            GM_LOGGER.error("error writing CSV output file: {}", e.toString());
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
                    GM_LOGGER.error("invalid items: {}", line);
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
                    germlineSVs.add(CsvFileWriter.fromString(line, true));
                }
                catch(Exception e)
                {
                    GM_LOGGER.error("sample({}) failed to parse line: {}", sampleId, line);
                }

                line = fileReader.readLine();
            }
        }
        catch (final IOException e)
        {
            GM_LOGGER.error("error reading file({}): {}", processedFile, e.toString());
        }
    }

    private void reprocessVariants(final String sampleId, final List<GermlineSV> germlineSVs)
    {
        /*
        GM_LOGGER.info("sample({}) reprocessing {} germline SVs for disruptions", sampleId, germlineSVs.size());

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
            GM_LOGGER.debug("sample({}) has {} disruptive SVs from total({})",
                    germlineSVs.get(0).SampleId, disruptiveSvCount, germlineSVs.size());
        }

        germlineSVs.forEach(x -> writeCsv(x));
        */
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

        setLogLevel(cmd);

        GermlineVcfReader germlineVcfReader = new GermlineVcfReader(cmd);
        germlineVcfReader.run();

        GM_LOGGER.info("VCF processing complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
