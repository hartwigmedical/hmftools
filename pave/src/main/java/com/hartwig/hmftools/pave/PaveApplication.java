package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.pave.PaveConfig.PON_ARTEFACTS_FILE;
import static com.hartwig.hmftools.pave.PaveConfig.PON_FILE;
import static com.hartwig.hmftools.pave.PaveConfig.PON_FILTERS;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveUtils.createRightAlignedVariant;
import static com.hartwig.hmftools.pave.annotation.PonAnnotation.PON_ARTEFACT_FILTER;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;

import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.VcfFileReader;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pave.annotation.Blacklistings;
import com.hartwig.hmftools.pave.annotation.ClinvarAnnotation;
import com.hartwig.hmftools.pave.annotation.GnomadAnnotation;
import com.hartwig.hmftools.pave.annotation.Mappability;
import com.hartwig.hmftools.pave.annotation.PonAnnotation;
import com.hartwig.hmftools.pave.annotation.Reportability;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class PaveApplication
{
    private final PaveConfig mConfig;
    private final ImpactClassifier mImpactClassifier;
    private final VariantImpactBuilder mImpactBuilder;
    private final GeneDataCache mGeneDataCache;
    private final GnomadAnnotation mGnomadAnnotation;

    private final PonAnnotation mPon;
    private final PonAnnotation mPonArtefacts;
    private final Mappability mMappability;
    private final ClinvarAnnotation mClinvar;
    private final Blacklistings mBlacklistings;
    private Reportability mReportability;

    private VcfWriter mVcfWriter;
    private final TranscriptWriter mTranscriptWriter;

    public PaveApplication(final CommandLine cmd)
    {
        final VersionInfo version = new VersionInfo("pave.version");
        PV_LOGGER.info("Pave version: {}", version.version());

        mConfig = new PaveConfig(cmd);

        mGeneDataCache = new GeneDataCache(
                cmd.getOptionValue(ENSEMBL_DATA_DIR), mConfig.RefGenVersion,
                cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION), true);

        mGnomadAnnotation = new GnomadAnnotation(cmd);

        mPon = new PonAnnotation(cmd.getOptionValue(PON_FILE), true);
        mPon.loadFilters(cmd.getOptionValue(PON_FILTERS));

        mPonArtefacts = new PonAnnotation(cmd.getOptionValue(PON_ARTEFACTS_FILE), false);

        mMappability = new Mappability(cmd);
        mClinvar = new ClinvarAnnotation(cmd);
        mBlacklistings = new Blacklistings(cmd);

        mImpactBuilder = new VariantImpactBuilder(mGeneDataCache);

        RefGenomeInterface refGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));

        mImpactClassifier = new ImpactClassifier(refGenome);

        mReportability = null;

        mVcfWriter = null;
        initialiseVcfWriter();

        mTranscriptWriter = new TranscriptWriter(mConfig);

        try
        {
            version.write(mConfig.OutputDir);
        }
        catch(IOException e)
        {
            System.exit(1);
        }
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            PV_LOGGER.error("invalid config, exiting");
            System.exit(1);
        }

        if(!mGeneDataCache.loadCache(mConfig.OnlyCanonical, false))
        {
            PV_LOGGER.error("gene data cache loading failed, exiting");
            System.exit(1);
        }

        mReportability = new Reportability(mGeneDataCache.getDriverPanel());

        if((mPon.isEnabled() && !mPon.hasValidData()) || (mPonArtefacts.isEnabled() && !mPonArtefacts.hasValidData()))
        {
            PV_LOGGER.error("invalid PON files, exiting");
            System.exit(1);
        }

        if(!mMappability.hasValidData())
        {
            PV_LOGGER.error("invalid mappability data, exiting");
            System.exit(1);
        }

        if(!mClinvar.hasValidData())
        {
            PV_LOGGER.error("invalid Clinvar data, exiting");
            System.exit(1);
        }

        if(!mBlacklistings.hasValidData())
        {
            PV_LOGGER.error("invalid blacklistings data, exiting");
            System.exit(1);
        }

        processVcfFile(mConfig.SampleId);

        mTranscriptWriter.close();

        PV_LOGGER.info("sample({}) annotation complete", mConfig.SampleId);
    }

    private void processVcfFile(final String sampleId)
    {
        PV_LOGGER.info("sample({}) reading VCF file({})", sampleId, mConfig.VcfFile);

        int variantCount = 0;

        VcfFileReader vcfFileReader = new VcfFileReader(mConfig.VcfFile);

        if(!vcfFileReader.fileValid())
        {
            PV_LOGGER.error("invalid somatic VCF file({})", mConfig.VcfFile);
            System.exit(1);
        }

        for(VariantContext variantContext : vcfFileReader.iterator())
        {
            if(!mConfig.SpecificRegions.isEmpty())
            {
                if(mConfig.SpecificRegions.stream().noneMatch(x -> x.containsPosition(variantContext.getContig(), variantContext.getStart())))
                    continue;
            }

            processVariant(variantContext);
            ++variantCount;

            if(variantCount > 0 && (variantCount % 100000) == 0)
            {
                PV_LOGGER.info("processed {} variants", variantCount);
            }
        }

        processPhasedVariants(NO_LOCAL_PHASE_SET);

        PV_LOGGER.info("sample({}) processed {} variants", sampleId, variantCount);

        mVcfWriter.close();
    }

    private void processVariant(final VariantContext variantContext)
    {
        if(!HumanChromosome.contains(variantContext.getContig()))
            return;

        VariantData variant = VariantData.fromContext(variantContext);

        if(mConfig.ReadPassOnly)
        {
            if(!variantContext.getFilters().isEmpty() && !variantContext.getFilters().contains(PASS_FILTER))
                return;
        }

        try
        {
            variant.setRealignedVariant(createRightAlignedVariant(variant, mImpactClassifier.refGenome()));

            findVariantImpacts(variant, mImpactClassifier, mGeneDataCache);

            processPhasedVariants(variant.localPhaseSet());

            if(!variant.hasLocalPhaseSet())
                processVariant(variant);
        }
        catch(Exception e)
        {
            PV_LOGGER.error("error processing var({})", variant);
            e.printStackTrace();
            System.exit(1);
        }
    }

    private void processPhasedVariants(int currentLocalPhaseSet)
    {
        List<VariantData> variants = mImpactClassifier.processPhasedVariants(currentLocalPhaseSet);

        if(variants != null)
            variants.forEach(x -> processVariant(x));
    }

    private void processVariant(final VariantData variant)
    {
        // can be null if no impacts exist for any transcript
        VariantImpact variantImpact = mImpactBuilder.createVariantImpact(variant);

        ponAnnotateAndFilter(variant);

        if(mConfig.SetReportable)
            mReportability.setReportability(variant, variantImpact);

        if(mConfig.WritePassOnly && !variant.filters().isEmpty())
            return;

        mVcfWriter.writeVariant(variant.context(), variant, variantImpact);

        if(!mConfig.WriteTranscriptFile)
            return;

        for(Map.Entry<String,List<VariantTransImpact>> entry : variant.getImpacts().entrySet())
        {
            final String geneName = entry.getKey();
            mTranscriptWriter.writeVariantData(variant, geneName);
        }
    }

    private void ponAnnotateAndFilter(final VariantData variant)
    {
        mGnomadAnnotation.annotateVariant(variant);
        mMappability.annotateVariant(variant);
        mClinvar.annotateVariant(variant);
        mBlacklistings.annotateVariant(variant);

        mPon.annotateVariant(variant);

        if(mPonArtefacts.getPonData(variant) != null)
            variant.addFilter(PON_ARTEFACT_FILTER);
    }

    public static void findVariantImpacts(
            final VariantData variant, final ImpactClassifier impactClassifier, final GeneDataCache geneDataCache)
    {
        boolean processed = false;

        List<GeneData> geneCandidates = geneDataCache.findGenes(variant.Chromosome, variant.Position, variant.EndPosition);

        if(!geneCandidates.isEmpty())
        {
            // analyse against each of the genes and their transcripts
            for(GeneData geneData : geneCandidates)
            {
                List<TranscriptData> transDataList =
                        geneDataCache.findTranscripts(geneData.GeneId, variant.Position, variant.EndPosition);

                // non-coding transcripts are skipped for now
                if(transDataList.isEmpty())
                    continue;

                for(TranscriptData transData : transDataList)
                {
                    VariantTransImpact transImpact = impactClassifier.classifyVariant(variant, transData);
                    processed = true;

                    // check right-alignment if the variant has microhomology
                    if(variant.realignedVariant() != null)
                    {
                        VariantTransImpact raTransImpact = impactClassifier.classifyVariant(variant.realignedVariant(), transData);

                        if(raTransImpact != null)
                        {
                            variant.realignedVariant().addImpact(geneData.GeneName, raTransImpact);
                            transImpact = ImpactClassifier.selectAlignedImpacts(transImpact, raTransImpact);
                        }
                    }

                    if(transImpact != null)
                        variant.addImpact(geneData.GeneName, transImpact);
                }
            }
        }

        // ensure all phased variants are cached
        if(!processed && variant.hasLocalPhaseSet())
            impactClassifier.phasedVariants().checkAddVariant(variant);
    }

    private void initialiseVcfWriter()
    {
        final VersionInfo version = new VersionInfo("pave.version");

        // append 'pave' to the input vcf file name if not specified
        String outputVcfFilename;

        if(mConfig.OutputVcfFile != null)
        {
            outputVcfFilename = mConfig.OutputVcfFile; // assumes includes path
        }
        else
        {
            String[] fileItems = mConfig.VcfFile.split("/");
            String filename = fileItems[fileItems.length - 1];
            int extensionIndex = filename.indexOf(".vcf");
            outputVcfFilename = mConfig.OutputDir + filename.substring(0, extensionIndex) + ".pave" + filename.substring(extensionIndex);

            if(!outputVcfFilename.endsWith(".gz")) // always writes zipped VCF even if input VCF isn't zipped
                outputVcfFilename += ".gz";
        }

        PV_LOGGER.info("writing VCF file({})", outputVcfFilename);

        mVcfWriter = new VcfWriter(outputVcfFilename, mConfig.VcfFile);

        mVcfWriter.writeHeader(
                version.version(), mGnomadAnnotation.hasData(), mPon.isEnabled(), mMappability.hasData(),
                mClinvar.hasData(), mBlacklistings.hasData(), mConfig.SetReportable);
    }


    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = PaveConfig.createOptions();

        try
        {
            final CommandLine cmd = createCommandLine(args, options);

            setLogLevel(cmd);

            PaveApplication paveApplication = new PaveApplication(cmd);
            paveApplication.run();
        }
        catch(ParseException e)
        {
            PV_LOGGER.warn(e);
            final HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("PaveApplication", options);
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
