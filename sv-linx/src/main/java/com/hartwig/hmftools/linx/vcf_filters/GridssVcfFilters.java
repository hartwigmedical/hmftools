package com.hartwig.hmftools.linx.vcf_filters;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_PAIR;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.AS;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.BEID;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.BEIDL;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineFilters.CAS;
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

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.filter.AlwaysPassFilter;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

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

public class GridssVcfFilters
{
    private final String mSampleId;
    private final String mOutputDir;
    private final String mScope;
    private final String mVcfFile;

    private BufferedWriter mCsvWriter;

    private final GermlineFilters mFilter;
    private final StructuralVariantFactory mSvFactory;
    private final FilterConfig mFilterConfig;
    private final boolean mLinkByAssembly;

    private SvGeneTranscriptCollection mGeneCollection;

    private List<AssemblyData> mSvAssemblyData;
    private Map<String,String> mSvFilterReasons;

    private static final String GRIDSS_VCF_FILE = "vcf";
    private static final String SCOPE = "scope";
    private static final String SAMPLE = "sample";
    private static final String REF_GENOME = "ref_genome";
    private static final String GENE_TRANSCRIPTS_DIR = "gene_trans_dir";
    private static final String LINK_BY_ASSEMBLY = "link_by_assembly";
    private static final String OUTPUT_DIR = "output_dir";
    private static final String LOG_DEBUG = "log_debug";

    private static final Logger LOGGER = LogManager.getLogger(GridssVcfFilters.class);

    public GridssVcfFilters(final CommandLine cmd)
    {
        mSampleId = cmd.getOptionValue(SAMPLE);
        mOutputDir = cmd.getOptionValue(OUTPUT_DIR);

        mVcfFile = cmd.getOptionValue(GRIDSS_VCF_FILE);
        mScope = cmd.getOptionValue(SCOPE);

        mFilterConfig = new FilterConfig(cmd);

        mFilter = new GermlineFilters(mFilterConfig);
        mSvFactory = new StructuralVariantFactory(new AlwaysPassFilter());

        mSvFilterReasons = Maps.newHashMap();

        mSvAssemblyData = Lists.newArrayList();
        mLinkByAssembly = cmd.hasOption(LINK_BY_ASSEMBLY);

        if(cmd.hasOption(GENE_TRANSCRIPTS_DIR))
        {
            mGeneCollection = new SvGeneTranscriptCollection();
            mGeneCollection.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));
            mGeneCollection.loadEnsemblData(true);
        }
        else
        {
            mGeneCollection = null;
        }
    }

    public void run()
    {
        LOGGER.info("scope({}) processing VCF({})", mScope, mVcfFile);

        processVcf();

        if(mLinkByAssembly)
            annotateAssembledLinks();

        writeSVs();
    }

    private void processVcf()
    {
        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                    mVcfFile, new VCFCodec(), false);

            reader.iterator().forEach(x -> processVariant(x));
        }
        catch(IOException e)
        {
            LOGGER.error("error reading vcf({}): {}", mVcfFile, e.toString());
        }
    }

    private void processVariant(final VariantContext variant)
    {
        LOGGER.trace("id({}) position({}: {})", variant.getID(), variant.getContig(), variant.getStart());

        // early exit
        if(mFilterConfig.RequirePass && !mFilterConfig.LogFiltered)
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

        final String filterStr = applyFilters(sv, variant);

        if(!filterStr.equals(PASS) && !mFilterConfig.LogFiltered)
        {
            purgeLastSv();
            return;
        }

        mSvFilterReasons.put(sv.id(), filterStr);

        if(mLinkByAssembly)
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
        if(!variant.getFilters().isEmpty() && mFilterConfig.RequirePass)
            return "GRIDSS_FILTERED";

        if(GermlineFilters.isImprecise(variant))
            return "IMPRECISE";

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

    private void annotationWithGenes(final StructuralVariant sv, final String[] geneAnnotations)
    {
        if(mGeneCollection == null)
            return;

        for(int be = SE_START; be <= SE_END; ++be)
        {
            if(sv.type() == SGL && be == SE_END)
                continue;

            final List<EnsemblGeneData> matchedGenes = mGeneCollection.findGenes(sv.chromosome(isStart(be)), sv.position(isStart(be)), 0);

            if(!matchedGenes.isEmpty())
            {
                String genesStr = "";

                for(final EnsemblGeneData gene : matchedGenes)
                {
                    genesStr = appendStr(genesStr, gene.GeneName, ';');
                }

                geneAnnotations[be] = genesStr;
            }
        }
    }

    private void annotateAssembledLinks()
    {
        for(int i = 0; i < mSvAssemblyData.size() - 1; ++i)
        {
            AssemblyData asmData1 = mSvAssemblyData.get(i);
            boolean linked = false;

            for(int j = i+1; j < mSvAssemblyData.size(); ++j)
            {
                AssemblyData asmData2 = mSvAssemblyData.get(j);

                for(int se1 = SE_START; se1 <= SE_END; ++se1)
                {
                    for(int se2 = SE_START; se2 <= SE_END; ++se2)
                    {
                        if(asmData1.hasMatch(asmData2, se1, se2))
                        {
                            asmData1.setLinkedData(asmData2.VcfId, isStart(se1));
                            asmData2.setLinkedData(asmData1.VcfId, isStart(se2));
                            linked = true;
                            break;
                        }
                    }

                    if(linked)
                        break;
                }

                if(linked)
                    break;
            }
        }
    }

    private void populateAssemblyLinks(final StructuralVariant sv, final String[] asmbLinks)
    {
        final AssemblyData asmData = mSvAssemblyData.stream().filter(x -> x.VcfId.equals(sv.id())).findFirst().orElse(null);

        if(asmData == null)
            return;

        asmbLinks[SE_START] = asmData.getLinkedSvIds()[SE_START];
        asmbLinks[SE_END] = asmData.getLinkedSvIds()[SE_END];
    }
    private void writeSVs()
    {
        if(mOutputDir.isEmpty())
            return;

        for(final StructuralVariant sv : mSvFactory.results())
        {
            final String filterStr = mSvFilterReasons.get(sv.id());
            writeCsv(sv, filterStr);
        }

        closeBufferedWriter(mCsvWriter);
    }

    private void writeCsv(final StructuralVariant sv, final String filter)
    {
        try
        {
            if(mCsvWriter == null)
            {
                String outputFileName = mOutputDir + mSampleId + "_filtered_germline_svs.csv";

                mCsvWriter = createBufferedWriter(outputFileName, false);

                mCsvWriter.write("SampleId,Id,Filter,GridssFilter,QualScore,Type");
                mCsvWriter.write(",ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd");
                // mCsvWriter.write(",NormalReadDepth,NormalAltCount,TumorReadDepth,TumorAltCount");
                mCsvWriter.write(",NormalREF,NormalRP,NormalRPQ,NormalSR,NormalSRQ,NormalVF");
                mCsvWriter.write(",GenesStart,GenesEnd,AsmbStart,AsmbEnd");

                mCsvWriter.newLine();
            }

            final VariantContext variant = sv.startContext();

            final Genotype normalGenotype = variant.getGenotype(0);
            // final Genotype tumorGenotype = variant.getGenotype(1);

            mCsvWriter.write(String.format("%s,%s,%s,%s,%f,%s",
                    mSampleId, sv.id(), filter, sv.filter(),
                    getDoubleValue(normalGenotype, QUAL), sv.type()));

            boolean sglBreakend = sv.type() == SGL;

            mCsvWriter.write(String.format(",%s,%s,%d,%d,%d,%d",
                    sv.chromosome(true), !sglBreakend ? sv.chromosome(false) : 0,
                    sv.position(true), !sglBreakend ? sv.position(false) : -1,
                    sv.orientation(true), !sglBreakend ? sv.orientation(false) : 0));

            // mCsvWriter.write(String.format(",%d,%d,%d,%d",
            //        normalGenotype.getDP(), normalGenotype.getAD(), tumorGenotype.getDP(), tumorGenotype.getAD()));

            mCsvWriter.write(String.format(",%d,%d,%.2f,%d,%.2f,%d",
                    getIntValue(normalGenotype, REF), getIntValue(normalGenotype, RP), getDoubleValue(normalGenotype, RPQ),
                    getIntValue(normalGenotype, SR), getDoubleValue(normalGenotype, SRQ), getIntValue(normalGenotype, VF)));

            String[] geneAnnotations = {"",""};
            annotationWithGenes(sv, geneAnnotations);

            String[] asmbSvIds = {"",""};
            populateAssemblyLinks(sv, asmbSvIds);

            mCsvWriter.write(String.format(",%s,%s,%s,%s",
                    geneAnnotations[SE_START], geneAnnotations[SE_END],
                    asmbSvIds[SE_START], asmbSvIds[SE_END]));

            mCsvWriter.newLine();
        }
        catch (final IOException e)
        {
            LOGGER.error("error writing CSV output file: {}", e.toString());
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
        final Options options = createBasicOptions();
        final CommandLine cmd = createCommandLine(args, options);

        if(cmd.hasOption(LOG_DEBUG))
            Configurator.setRootLevel(Level.DEBUG);

        GridssVcfFilters gridssVcfFilters = new GridssVcfFilters(cmd);

        gridssVcfFilters.run();

        LOGGER.info("VCF processing complete");
    }

    @NotNull
    private static Options createBasicOptions()
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the tumor sample");
        options.addOption(GRIDSS_VCF_FILE, true, "Path to the GRIDSS structural variant VCF file");
        options.addOption(REF_GENOME, true, "Ref genome fasta file");
        options.addOption(GENE_TRANSCRIPTS_DIR, true, "Ensembl data cache directory");
        options.addOption(LINK_BY_ASSEMBLY, false, "Look for assembled links");
        options.addOption(SCOPE, true, "Scope: germline or somatic");
        options.addOption(OUTPUT_DIR, true, "Path to write results");
        options.addOption(LOG_DEBUG, false, "Log verbose");

        FilterConfig.addCommandLineOptions(options);

        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
