package com.hartwig.hmftools.pave;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.VariantConsequence.VARIANT_CONSEQ_DELIM;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.PaveConstants.DELIM;
import static com.hartwig.hmftools.pave.PaveUtils.createRightAlignedVariant;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;
import static com.hartwig.hmftools.pave.compare.ComparisonUtils.effectsMatch;
import static com.hartwig.hmftools.pave.compare.ComparisonUtils.ignoreSnpEffAnnotation;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.version.VersionInfo;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationParser;
import com.hartwig.hmftools.pave.compare.ComparisonUtils;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;

public class PaveApplication
{
    private final PaveConfig mConfig;
    private final ImpactClassifier mImpactClassifier;
    private final VariantImpactBuilder mImpactBuilder;
    private final GeneDataCache mGeneDataCache;

    // for comparison with SnpEff
    private int mTransEffectMatched;
    private int mTransTotalComparisons;

    private VcfWriter mVcfWriter;
    private BufferedWriter mCsvTranscriptWriter;

    public PaveApplication(final CommandLine cmd)
    {
        mConfig = new PaveConfig(cmd);

        mGeneDataCache = new GeneDataCache(
                cmd.getOptionValue(ENSEMBL_DATA_DIR), mConfig.RefGenVersion,
                cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION), mConfig.CompareSnpEff, true);

        mImpactBuilder = new VariantImpactBuilder(mGeneDataCache);

        RefGenomeInterface refGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));

        mImpactClassifier = new ImpactClassifier(refGenome);

        mVcfWriter = null;
        initialiseVcfWriter();

        if(mConfig.WriteTranscriptCsv)
            initialiseTranscriptWriter();

        mTransEffectMatched = 0;
        mTransTotalComparisons = 0;
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            PV_LOGGER.error("invalid config, exiting");
            System.exit(1);
        }

        if(!mGeneDataCache.loadCache())
        {
            PV_LOGGER.error("Ensembl data cache loading failed, exiting");
            System.exit(1);
        }

        processVcfFile(mConfig.SampleId);

        closeBufferedWriter(mCsvTranscriptWriter);

        PV_LOGGER.info("sample({}) annotation complete", mConfig.SampleId);
    }

    private void processVcfFile(final String sampleId)
    {
        PV_LOGGER.info("sample({}) reading VCF file({})", sampleId, mConfig.VcfFile);

        int variantCount = 0;

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                    mConfig.VcfFile, new VCFCodec(), false);

            for(VariantContext variantContext : reader.iterator())
            {
                processVariant(variantContext);
                ++variantCount;

                if(variantCount > 0 && (variantCount % 10000) == 0)
                {
                    PV_LOGGER.debug("processed {} variants", variantCount);
                }
            }

            processPhasedVariants(NO_LOCAL_PHASE_SET);
        }
        catch(IOException e)
        {
            PV_LOGGER.error(" failed to read somatic VCF file({}): {}", mConfig.VcfFile, e.toString());
        }

        PV_LOGGER.info("sample({}) processed {} variants", sampleId, variantCount);

        if(mConfig.CompareSnpEff)
        {
            PV_LOGGER.info("sample({}) transcript comparison total({} matched={} diffs={})",
                    sampleId, mTransTotalComparisons, mTransEffectMatched, mTransTotalComparisons - mTransEffectMatched);
        }

        mVcfWriter.close();
    }

    private void processVariant(final VariantContext variantContext)
    {
        VariantData variant = VariantData.fromContext(variantContext);

        variant.setRealignedVariant(createRightAlignedVariant(variant, mImpactClassifier.refGenome()));

        findVariantImpacts(variant, mImpactClassifier, mGeneDataCache);

        processPhasedVariants(variant.localPhaseSet());

        if(!variant.hasLocalPhaseSet())
            processVariant(variant);
    }

    private void processPhasedVariants(int currentLocalPhaseSet)
    {
        List<VariantData> variants = mImpactClassifier.processPhasedVariants(currentLocalPhaseSet);

        if(variants != null)
            variants.forEach(x -> processVariant(x));
    }

    private void processVariant(final VariantData variant)
    {
        VariantImpact variantImpact = mImpactBuilder.createVariantImpact(variant);

        mVcfWriter.writeVariant(variant.context(), variant, variantImpact);

        if(!mConfig.WriteTranscriptCsv)
            return;

        if(!mConfig.CompareSnpEff)
        {
            for(Map.Entry<String,List<VariantTransImpact>> entry : variant.getImpacts().entrySet())
            {
                final String geneName = entry.getKey();
                writeVariantTranscriptData(variant, geneName);
            }
        }
        else
        {
            checkAndWriteTransDifferences(variant);
        }
    }

    public static void findVariantImpacts(
            final VariantData variant, final ImpactClassifier impactClassifier, final GeneDataCache geneDataCache)
    {
        List<GeneData> geneCandidates = geneDataCache.findGenes(variant.Chromosome, variant.Position, variant.EndPosition);

        if(geneCandidates.isEmpty())
            return;

        // analyse against each of the genes and their transcripts
        for(GeneData geneData : geneCandidates)
        {
            List<TranscriptData> transDataList = geneDataCache.findTranscripts(geneData.GeneId, variant.Position, variant.EndPosition);

            // non-coding transcripts are skipped for now
            if(transDataList.isEmpty())
                continue;

            for(TranscriptData transData : transDataList)
            {
                VariantTransImpact transImpact = impactClassifier.classifyVariant(variant, transData);

                // check right-alignment if the variant has microhomology
                if(variant.realignedVariant() != null)
                {
                    VariantTransImpact raTransImpact = impactClassifier.classifyVariant(variant.realignedVariant(), transData);
                    transImpact = ImpactClassifier.selectAlignedImpacts(transImpact, raTransImpact);
                }

                if(transImpact != null)
                    variant.addImpact(geneData.GeneName, transImpact);
            }
        }
    }

    private void initialiseVcfWriter()
    {
        final VersionInfo version = new VersionInfo("pave.version");

        String vcfFilename = mConfig.OverwriteVcf ?
                mConfig.VcfFile : mConfig.OutputDir + mConfig.SampleId + ".pave.annotated.vcf";

        mVcfWriter = new VcfWriter(vcfFilename, mConfig.VcfFile);

        // TO-DO
        // mVcfWriter.writeHeader(version.version());
        mVcfWriter.writeHeader("1.0");
    }

    private void initialiseTranscriptWriter()
    {
        try
        {
            String fileSuffix = mConfig.CompareSnpEff ? ".pave.transcript_compare.csv" : ".pave.transcript.csv";
            String transFileName = mConfig.OutputDir + mConfig.SampleId + fileSuffix;
            mCsvTranscriptWriter = createBufferedWriter(transFileName, false);

            StringJoiner sj = new StringJoiner(DELIM);
            sj.add(VariantData.csvCommonHeader());
            sj.add("GeneId,GeneName");
            sj.add(VariantTransImpact.csvHeader());
            sj.add(CodingContext.csvHeader());
            sj.add(ProteinContext.csvHeader());

            if(mConfig.CompareSnpEff)
                sj.add("SnpEffEffects,SnpEffHgvsCoding,SnpEffHgvsProtein");

            mCsvTranscriptWriter.write(sj.toString());
            mCsvTranscriptWriter.newLine();
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return;
        }
    }

    private void writeVariantTranscriptData(final VariantData variant, final String geneName)
    {
        if(mCsvTranscriptWriter == null)
            return;

        List<VariantTransImpact> geneImpacts = variant.getImpacts().get(geneName);

        if(geneImpacts == null)
            return;

        try
        {
            for(VariantTransImpact impact : geneImpacts)
            {
                if(impact.TransData == null)
                    continue;

                mCsvTranscriptWriter.write(String.format("%s", variant.toCsv()));
                mCsvTranscriptWriter.write(String.format(",%s,%s,%s,%s,%s",
                        impact.TransData.GeneId, geneName, impact.toCsv(), impact.codingContext().toCsv(),
                        impact.proteinContext() != null ? impact.proteinContext().toCsv() : ProteinContext.empty()));
                mCsvTranscriptWriter.newLine();
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to write variant CSV file: {}", e.toString());
            return;
        }
    }

    private void checkAndWriteTransDifferences(final VariantData variant)
    {
        // extract SnpEff data for comparison sake
        final VariantContext variantContext = variant.context();

        List<SnpEffAnnotation> snpEffAnnotations = SnpEffAnnotationParser.fromContext(variantContext).stream()
                .filter(x -> !ignoreSnpEffAnnotation(x))
                .collect(Collectors.toList());

        if(variant.getImpacts().isEmpty())
        {
            if(snpEffAnnotations != null && !snpEffAnnotations.isEmpty())
            {
                for(SnpEffAnnotation annotation : snpEffAnnotations)
                {
                    GeneData geneData = mGeneDataCache.findSnpEffGeneData(annotation.geneID(), annotation.gene());

                    if(geneData == null)
                    {
                        //  PV_LOGGER.debug("ignoring unknown gene({}:{})", annotation.geneID(), annotation.gene());
                        continue;
                    }

                    writeVariantTranscriptData(variant, geneData.GeneName, snpEffAnnotations);
                }
            }
            else
            {
                writeVariantTranscriptData(variant, null, snpEffAnnotations);
            }

            return;
        }

        // analyse against each of the genes and their transcripts
        for(Map.Entry<String,List<VariantTransImpact>> entry : variant.getImpacts().entrySet())
        {
            final String geneName = entry.getKey();
            final String snpEffGeneName = mGeneDataCache.getSnpEffGeneName(geneName);

            writeVariantTranscriptData(
                    variant, geneName,
                    snpEffAnnotations.stream().filter(x -> x.gene().equals(snpEffGeneName)).collect(Collectors.toList()));
        }
    }

    private void writeVariantTranscriptData(final VariantData variant, final String geneName, final List<SnpEffAnnotation> annotations)
    {
        if(mCsvTranscriptWriter == null)
            return;

        try
        {
            List<String> transcriptLines = Lists.newArrayList();

            List<SnpEffAnnotation> processedAnnotations = Lists.newArrayList();

            List<VariantTransImpact> geneImpacts = variant.getImpacts().get(geneName);

            if(geneImpacts != null)
            {
                for(VariantTransImpact impact : geneImpacts)
                {
                    if(impact.TransData == null)
                        continue;

                    List<SnpEffAnnotation> matchedAnnotations = ComparisonUtils.findMatchingAnnotations(impact, annotations);
                    processedAnnotations.addAll(matchedAnnotations);

                    ++mTransTotalComparisons;

                    if(!matchedAnnotations.isEmpty() && matchedAnnotations.stream().anyMatch(x -> effectsMatch(impact, x)))
                    {
                        ++mTransEffectMatched;

                        if(mConfig.WriteDiffs)
                            continue;
                    }

                    StringJoiner sj = new StringJoiner(DELIM);
                    sj.add(variant.toCsv());

                    sj.add(impact.TransData.GeneId);
                    sj.add(geneName);


                    sj.add(impact.toCsv());
                    sj.add(impact.codingContext().toCsv());
                    sj.add(impact.proteinContext() != null ? impact.proteinContext().toCsv() : ProteinContext.empty());

                    if(!matchedAnnotations.isEmpty())
                    {
                        if(matchedAnnotations.size() == 1)
                        {
                            sj.add(matchedAnnotations.get(0).effects());
                            sj.add(matchedAnnotations.get(0).hgvsCoding());
                            sj.add(matchedAnnotations.get(0).hgvsProtein());
                        }
                        else
                        {
                            StringJoiner maSj = new StringJoiner(VARIANT_CONSEQ_DELIM);
                            matchedAnnotations.forEach(x -> maSj.add(x.effects()));
                            sj.add(maSj.toString());

                            StringJoiner maC = new StringJoiner(VARIANT_CONSEQ_DELIM);
                            matchedAnnotations.stream().filter(x -> !x.hgvsCoding().isEmpty()).forEach(x -> maC.add(x.hgvsCoding()));
                            sj.add(maC.toString());

                            StringJoiner maP = new StringJoiner(VARIANT_CONSEQ_DELIM);
                            matchedAnnotations.stream().filter(x -> !x.hgvsProtein().isEmpty()).forEach(x -> maP.add(x.hgvsProtein()));
                            sj.add(maP.toString());
                        }
                    }
                    else
                    {
                        sj.add("UNMATCHED,na,na");
                    }

                    transcriptLines.add(sj.toString());
                }
            }

            for(SnpEffAnnotation annotation : annotations)
            {
                if(processedAnnotations.contains(annotation))
                    continue;

                final TranscriptData transData = mGeneDataCache.getTranscriptData(annotation.geneID(), annotation.featureID());

                if(transData == null)
                    continue;

                ++mTransTotalComparisons;

                StringJoiner sj = new StringJoiner(DELIM);
                sj.add(variant.toCsv());
                sj.add(annotation.geneID());
                sj.add(annotation.gene());
                sj.add(annotation.featureID());
                sj.add("UNMATCHED");
                sj.add(annotation.effects());
                transcriptLines.add(sj.toString());
            }

            for(String transData : transcriptLines)
            {
                mCsvTranscriptWriter.write(transData);
                mCsvTranscriptWriter.newLine();
            }
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to write variant CSV file: {}", e.toString());
            return;
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = PaveConfig.createOptions();
        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        PaveApplication paveApplication = new PaveApplication(cmd);
        paveApplication.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
