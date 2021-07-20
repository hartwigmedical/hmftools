package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.impact.ImpactAnnotator.findVariantImpacts;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.patientdb.dao.SomaticVariantDAO;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;
import org.jooq.Record;
import org.jooq.Result;

public class ImpactComparisons
{
    private final ImpactConfig mConfig;
    private final ImpactClassifier mImpactClassifier;
    private final VariantImpactBuilder mImpactBuilder;
    private final GeneDataCache mGeneDataCache;

    private BufferedWriter mCsvWriter;

    // counters
    private int mVariantCount;
    private int mTransVariantCount;
    private final PerformanceCounter mPerfCounter;

    public ImpactComparisons(final CommandLine cmd)
    {
        mConfig = new ImpactConfig(cmd);

        mGeneDataCache = new GeneDataCache(
                cmd.getOptionValue(ENSEMBL_DATA_DIR), mConfig.RefGenVersion, cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION));

        mImpactBuilder = new VariantImpactBuilder(mGeneDataCache);

        RefGenomeInterface refGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));
        mImpactClassifier = new ImpactClassifier(refGenome);

        mCsvWriter = null;

        mPerfCounter = new PerformanceCounter("ClassifyVariant");
        mVariantCount = 0;
        mTransVariantCount = 0;
    }

    public void close()
    {
        closeBufferedWriter(mCsvWriter);
    }

    public void run()
    {
        if(!mConfig.multiSampleValid())
        {
            SG_LOGGER.error("invalid config, exiting");
            System.exit(1);
        }

        if(!mGeneDataCache.loadCache())
        {
            SG_LOGGER.error("Ensembl data cache loading failed, exiting");
            System.exit(1);
        }

        initialiseImpactWriter();

        for(int i = 0; i < mConfig.SampleIds.size(); ++i)
        {
            String sampleId = mConfig.SampleIds.get(i);

            mPerfCounter.start();
            processSample(sampleId);
            mPerfCounter.stop();

            if(i > 0 && (i % 100) == 0)
            {
                SG_LOGGER.info("processing {} samples", i);
            }
        }

        close();

        mPerfCounter.logStats();

        SG_LOGGER.info("impact comparison complete");
    }

    private void processSample(final String sampleId)
    {
        mVariantCount = 0;
        mTransVariantCount = 0;

        Result<Record> result = mConfig.DbAccess.context().select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(SOMATICVARIANT.GENE.isNotNull())
                .orderBy(SOMATICVARIANT.CHROMOSOME,SOMATICVARIANT.POSITION)
                .fetch();

        for (Record record : result)
        {
            processVariant(sampleId, SomaticVariantDAO.buildFromRecord(record));
        }

        SG_LOGGER.debug("sample({}) processed {} variants and transcripts({})",
                sampleId, mVariantCount, mTransVariantCount);
    }

    private void processVariant(final String sampleId, final SomaticVariant somaticVariant)
    {
        ++mVariantCount;

        // generate variant impact data and then write comparison results to CSV file
        VariantData variant = new VariantData(
                somaticVariant.chromosome(), (int)somaticVariant.position(), somaticVariant.ref(), somaticVariant.alt());

        boolean phasedInframeIndel = somaticVariant.phasedInframeIndelIdentifier() != null && somaticVariant.phasedInframeIndelIdentifier() > 0;
        variant.setVariantDetails(phasedInframeIndel, somaticVariant.microhomology(), somaticVariant.repeatSequence());

        findVariantImpacts(variant, mImpactClassifier, mGeneDataCache);

        VariantImpact variantImpact = mImpactBuilder.createVariantImpact(variant);

        if(SG_LOGGER.isDebugEnabled())
        {
            if(variantImpact.CanonicalCodingEffect != somaticVariant.canonicalCodingEffect())
            {
                //SG_LOGGER.debug("sample({}) var({}) diff canonical coding: new({}) snpEff({})",
                //        sampleId, variant.toString(), variantImpact.CanonicalCodingEffect, somaticVariant.canonicalCodingEffect());
            }
        }

        writeVariantData(sampleId, variant, variantImpact, somaticVariant);
    }

    private void initialiseImpactWriter()
    {
        try
        {
            String fileName = mConfig.OutputDir + "VARIANT_IMPACT_COMPARE.csv";
            mCsvWriter = createBufferedWriter(fileName, false);

            mCsvWriter.write("SampleId,");
            mCsvWriter.write(VariantData.csvCommonHeader());
            mCsvWriter.write(",GeneId,GeneName,IsDriver,CanonEffect,CanonCodingEffect");
            mCsvWriter.write(",WorstEffect,WorstCodingEffect,WorstTrans,GenesAffected");
            mCsvWriter.write(",SnpEffGeneName,SnpEffCanonEffect,SnpEffCanonCodingEffect");
            mCsvWriter.write(",SnpEffWorstEffect,SnpEffWorstCodingEffect,SnpEffWorstTrans,SnpEffGenesAffected");
            mCsvWriter.newLine();
        }
        catch(IOException e)
        {
            SG_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return;
        }
    }

    private void writeVariantData(
            final String sampleId, final VariantData variant, final VariantImpact variantImpact, final SomaticVariant somaticVariant)
    {
        try
        {
            mCsvWriter.write(String.format("%s,%s,%s,%s",
                    sampleId, variant.toCsv(), variantImpact.CanonicalGeneId, variantImpact.CanonicalGeneName));

            boolean isDriver = mGeneDataCache.getDriverPanelGenes().contains(variantImpact.CanonicalGeneName);

            if(!isDriver && !variantImpact.CanonicalGeneName.equals(somaticVariant.gene()))
                isDriver = mGeneDataCache.getDriverPanelGenes().contains(somaticVariant.gene());

            mCsvWriter.write(String.format(",%s,%s,%s,%s,%s,%s,%d",
                    isDriver, variantImpact.CanonicalEffect, variantImpact.CanonicalCodingEffect,
                    variantImpact.WorstEffect, variantImpact.WorstCodingEffect, variantImpact.WorstTranscript, variantImpact.GenesAffected));

            mCsvWriter.write(String.format(",%s,%s,%s,%s,%s,%s,%d",
                    somaticVariant.gene(), somaticVariant.canonicalEffect(), somaticVariant.canonicalCodingEffect(), somaticVariant.worstEffect(),
                    somaticVariant.worstCodingEffect(), somaticVariant.worstEffectTranscript(), somaticVariant.genesAffected()));

            mCsvWriter.newLine();
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

        ImpactComparisons impactComparison = new ImpactComparisons(cmd);
        impactComparison.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
