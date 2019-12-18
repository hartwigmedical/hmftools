package com.hartwig.hmftools.linx.vcf_filters;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.BND;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.INV;
import static com.hartwig.hmftools.common.variant.structural.StructuralVariantType.SGL;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStr;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.linx.types.SvVarData.isStart;
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
import static com.hartwig.hmftools.linx.vcf_filters.GermlineVcfConfig.GENE_PANEL_FILE;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineVcfConfig.GENE_TRANSCRIPTS_DIR;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineVcfConfig.LOG_DEBUG;
import static com.hartwig.hmftools.linx.vcf_filters.GermlineVcfConfig.loadVcfFiles;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.FileVisitOption;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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

public class GermlineVcfReader
{
    private BufferedWriter mCsvWriter;

    private final GermlineFilters mFilter;
    private StructuralVariantFactory mSvFactory;
    private final GermlineVcfConfig mConfig;

    private final SvGeneTranscriptCollection mGeneCollection;
    private final Map<String,List<EnsemblGeneData>> mGenePanel;

    private final List<String> mVcfFiles;

    private List<AssemblyData> mSvAssemblyData;
    private Map<String,String> mSvFilterReasons;

    private static final Logger LOGGER = LogManager.getLogger(GermlineVcfReader.class);

    public GermlineVcfReader(final CommandLine cmd)
    {
        mConfig = new GermlineVcfConfig(cmd);
        mFilter = new GermlineFilters(mConfig);
        mSvFactory = null;

        mVcfFiles = Lists.newArrayList();
        mSvFilterReasons = Maps.newHashMap();
        mSvAssemblyData = Lists.newArrayList();

        mGenePanel = Maps.newHashMap();

        if(cmd.hasOption(GENE_TRANSCRIPTS_DIR))
        {
            mGeneCollection = new SvGeneTranscriptCollection();
            mGeneCollection.setDataPath(cmd.getOptionValue(GENE_TRANSCRIPTS_DIR));
            mGeneCollection.loadEnsemblData(true);

            if(cmd.hasOption(GENE_PANEL_FILE))
                loadGenePanel(cmd.getOptionValue(GENE_PANEL_FILE));
        }
        else
        {
            mGeneCollection = null;
        }
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

    private void loadGenePanel(final String geneFile)
    {
        // gene_panel.csv
        if (!Files.exists(Paths.get(geneFile)))
            return;

        try
        {
            List<String> fileContents = Files.readAllLines(new File(geneFile).toPath());

            for(final String gene : fileContents)
            {
                EnsemblGeneData geneData = mGeneCollection.getGeneDataByName(gene);

                if(geneData == null)
                {
                    LOGGER.warn("gene({}) not found in Ensembl data cache", gene);
                    continue;
                }

                List<EnsemblGeneData> chrGenesList = mGenePanel.get(geneData.Chromosome);

                if(chrGenesList == null)
                {
                    chrGenesList = Lists.newArrayList();
                    mGenePanel.put(geneData.Chromosome, chrGenesList);
                }

                chrGenesList.add(geneData);
            }

            LOGGER.info("loaded genePanelFile({}) with {} genes", geneFile, fileContents.size());
        }
        catch(IOException e)
        {
            LOGGER.error("failed to load gene panel file({}): {}", geneFile, e.toString());
        }
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
                annotateAssembledLinks();

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

        if(mConfig.RequireGene && !hasGeneAnnotation(sv))
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

    private String annotateWithGenePanel(final StructuralVariant sv)
    {
        if(sv.type() != DEL && sv.type() != DUP)
            return "";

        List<EnsemblGeneData> genesList = mGenePanel.get(sv.chromosome(true));

        if(genesList == null)
            return "";

        String genesStr = "";

        for(final EnsemblGeneData geneData : genesList)
        {
            if(sv.position(true) < geneData.GeneStart && sv.position(false) > geneData.GeneEnd)
            {
                genesStr = appendStr(genesStr, geneData.GeneName, ';');
            }
        }

        return genesStr;
    }

    private boolean hasGeneAnnotation(final StructuralVariant sv)
    {
        List<EnsemblGeneData> genesList = mGenePanel.get(sv.chromosome(true));

        if(genesList != null)
        {
            for (final EnsemblGeneData geneData : genesList)
            {
                // fully overlapping DEL or DUP
                if (sv.type() == DEL || sv.type() == DUP)
                {
                    if (sv.position(true) < geneData.GeneStart && sv.position(false) > geneData.GeneEnd)
                        return true;
                }

                // breakend falling in the gene
                if (sv.position(true) > geneData.GeneStart && sv.position(true) < geneData.GeneEnd)
                    return true;

                if (sv.type() == DEL || sv.type() == DUP | sv.type() == INV)
                {
                    if (sv.position(false) > geneData.GeneStart && sv.position(false) < geneData.GeneEnd)
                        return true;
                }
            }
        }

        // check other end of BND
        if(sv.type() == BND)
        {
            genesList = mGenePanel.get(sv.chromosome(false));

            if(genesList == null)
                return false;

            for (final EnsemblGeneData geneData : genesList)
            {
                if (sv.position(false) > geneData.GeneStart && sv.position(false) < geneData.GeneEnd)
                    return true;
            }
        }

        return false;
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
        for(final StructuralVariant sv : mSvFactory.results())
        {
            final String filterStr = mSvFilterReasons.get(sv.id());
            writeCsv(sv, filterStr);
        }
    }

    private void writeCsv(final StructuralVariant sv, final String filter)
    {
        try
        {
            if(mCsvWriter == null)
            {
                String outputFileName = mConfig.OutputDir;

                if(!mConfig.VcfFile.isEmpty())
                    outputFileName += mConfig.SampleId + "_filtered_germline_svs.csv";
                else
                    outputFileName += "LNX_GERMLINE_SVS.csv";

                mCsvWriter = createBufferedWriter(outputFileName, false);

                mCsvWriter.write("SampleId,Id,Filter,GridssFilter,QualScore,Type");
                mCsvWriter.write(",ChrStart,ChrEnd,PosStart,PosEnd,OrientStart,OrientEnd");
                // mCsvWriter.write(",NormalReadDepth,NormalAltCount,TumorReadDepth,TumorAltCount");
                mCsvWriter.write(",NormalREF,NormalRP,NormalRPQ,NormalSR,NormalSRQ,NormalVF");
                mCsvWriter.write(",InsertSequence,Homology");
                mCsvWriter.write(",GenesStart,GenesEnd,GenePanelOverlaps,AsmbStart,AsmbEnd");

                mCsvWriter.newLine();
            }

            final VariantContext variant = sv.startContext();

            final Genotype normalGenotype = variant.getGenotype(0);
            final String sampleName = normalGenotype.getSampleName();
            // final Genotype tumorGenotype = variant.getGenotype(1);

            mCsvWriter.write(String.format("%s,%s,%s,%s,%f,%s",
                    sampleName, sv.id(), filter, sv.filter(),
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

            final String homology = variant.getAttributeAsString(HOMSEQ, "");
            mCsvWriter.write(String.format(",%s,%s", sv.insertSequence(), homology));

            String[] geneAnnotations = {"",""};
            annotationWithGenes(sv, geneAnnotations);

            String[] asmbSvIds = {"",""};
            populateAssemblyLinks(sv, asmbSvIds);

            String genePanelOverlaps = annotateWithGenePanel(sv);

            mCsvWriter.write(String.format(",%s,%s,%s,%s,%s",
                    geneAnnotations[SE_START], geneAnnotations[SE_END], genePanelOverlaps,
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
