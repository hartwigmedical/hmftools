package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.VariantConsequence.INTRAGENIC_VARIANT;
import static com.hartwig.hmftools.common.variant.VariantConsequence.NON_CODING_TRANSCRIPT_VARIANT;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.impact.ImpactConfig.REF_GENOME;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotation;
import com.hartwig.hmftools.common.variant.snpeff.SnpEffAnnotationParser;

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

public class ImpactAnnotator
{
    private final ImpactConfig mConfig;
    private final ImpactClassifier mImpactClassifier;
    private final GeneDataCache mGeneDataCache;
    private BufferedWriter mCsvWriter;
    private BufferedWriter mCsvTranscriptWriter;

    private int mVariantCount;
    private int mVariantTransCount;
    private final PerformanceCounter mPerfCounter;

    public ImpactAnnotator(final CommandLine cmd)
    {
        mConfig = new ImpactConfig(cmd);

        mGeneDataCache = new GeneDataCache(cmd.getOptionValue(ENSEMBL_DATA_DIR), mConfig.RefGenVersion);

        RefGenomeInterface refGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));

        mImpactClassifier = new ImpactClassifier(refGenome);

        mPerfCounter = new PerformanceCounter("ClassifyVariant");
        mVariantCount = 0;
        mVariantTransCount = 0;

        mCsvWriter = null;
        mCsvTranscriptWriter = null;
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            SG_LOGGER.error("invalid config, exiting");
            System.exit(1);
        }

        if(!mGeneDataCache.loadCache())
        {
            SG_LOGGER.error("Ensembl data cache loading failed, exiting");
            System.exit(1);
        }

        if(mConfig.WriteCsv || mConfig.WriteTranscriptCsv)
            initialiseCsvWriters();

        SG_LOGGER.info("annotating variants with gene impacts");

        processVcfFile();

        closeBufferedWriter(mCsvWriter);
        closeBufferedWriter(mCsvTranscriptWriter);

        SG_LOGGER.info("processed {} variants and {} variant-transcript classifications", mVariantCount, mVariantTransCount);
        mPerfCounter.logStats();

        SG_LOGGER.info("annotation complete");
    }

    private void processVcfFile()
    {
        // CompoundFilter filter = new CompoundFilter(true);

        SG_LOGGER.info("sample({}) reading VCF file({})", mConfig.SampleId, mConfig.VcfFile);

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                    mConfig.VcfFile, new VCFCodec(), false);

            for (VariantContext variantContext : reader.iterator())
            {
                // if (!filter.test(variantContext))
                //    continue;

                mPerfCounter.start();
                ++mVariantCount;

                processVariant(variantContext);

                mPerfCounter.stop();
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error(" failed to read somatic VCF file({}): {}", mConfig.VcfFile, e.toString());
        }
    }

    private void processVariant(final VariantContext variantContext)
    {
        // SomaticVariantFactory variantFactory = new SomaticVariantFactory(filter);
        // final SomaticVariant variant = variantFactory.createVariant(mConfig.SampleId, variantContext).orElse(null);

        VariantData variant = VariantData.fromContext(variantContext);

        // extract SnpEff data for comparison sake
        List<SnpEffAnnotation> snpEffAnnotations = SnpEffAnnotationParser.fromContext(variantContext);
        variant.setSnpEffAnnotations(snpEffAnnotations);

        List<GeneData> geneCandidates = mGeneDataCache.findGenes(variant.Chromosome, variant.Position);

        if(geneCandidates.isEmpty())
        {
            // could skip these if sure that valid genes aren't being incorrectly missed
            variant.addImpact(new VariantTransImpact(null, INTRAGENIC_VARIANT));

            if(!snpEffAnnotations.isEmpty())
            {
                for(SnpEffAnnotation annotation : snpEffAnnotations)
                {
                    String geneId = annotation.geneID();

                    if(geneId.isEmpty())
                        continue;

                    GeneData geneData = mGeneDataCache.getEnsemblCache().getGeneDataById(geneId);

                    if(geneData == null)
                    {
                        SG_LOGGER.debug("ignoring unknown gene({}:{})", annotation.geneID(), annotation.gene());
                    }
                    else
                    {
                        writeVariantCsvData(variant, geneData);
                    }
                }
            }
            else
            {
                writeVariantCsvData(variant, null);
            }

            return;
        }

        // analyse against each of the genes and their transcripts
        for(GeneData geneData : geneCandidates)
        {
            List<TranscriptData> transDataList = mGeneDataCache.findTranscripts(geneData.GeneId, variant.Position);

            if(transDataList.isEmpty())
            {
                variant.addImpact(new VariantTransImpact(null, NON_CODING_TRANSCRIPT_VARIANT));
                continue;
            }

            for(TranscriptData transData : transDataList)
            {
                VariantTransImpact transImpact = mImpactClassifier.classifyVariant(variant, transData);

                if(transImpact != null)
                    variant.addImpact(transImpact);

                ++mVariantTransCount;
            }

            writeVariantCsvData(variant, geneData);
        }
    }

    private void initialiseCsvWriters()
    {
        try
        {
            if(mConfig.WriteCsv)
            {
                String fileName = mConfig.OutputDir + mConfig.SampleId + ".sage.variants.csv";
                mCsvWriter = createBufferedWriter(fileName, false);
                mCsvWriter.write(VariantData.csvHeader());
                mCsvWriter.newLine();
            }

            if(mConfig.WriteTranscriptCsv)
            {
                String transFileName = mConfig.OutputDir + mConfig.SampleId + ".sage.variant_transcripts.csv";
                mCsvTranscriptWriter = createBufferedWriter(transFileName, false);
                mCsvTranscriptWriter.write(VariantData.csvTranscriptHeader());
                mCsvTranscriptWriter.newLine();
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return;
        }
    }

    private void writeVariantCsvData(final VariantData variant, final GeneData geneData)
    {
        try
        {
            if(mCsvWriter != null)
            {
                mCsvWriter.write(variant.csvData(geneData));
                mCsvWriter.newLine();
            }

            if(mCsvTranscriptWriter != null)
            {
                for(String transData : variant.csvTranscriptData(geneData))
                {
                    mCsvTranscriptWriter.write(transData);
                    mCsvTranscriptWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to write variant CSV file: {}", e.toString());
            return;
        }

    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = ImpactConfig.createOptions();
        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        ImpactAnnotator impactAnnotator = new ImpactAnnotator(cmd);
        impactAnnotator.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
