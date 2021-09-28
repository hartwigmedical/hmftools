package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadSampleIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;
import static com.hartwig.hmftools.pave.PaveApplication.findVariantImpacts;
import static com.hartwig.hmftools.pave.PaveConfig.VI_LOGGER;
import static com.hartwig.hmftools.pave.compare.RefVariantData.hasCodingEffectDiff;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.PerformanceCounter;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.SomaticVariantDAO;
import com.hartwig.hmftools.pave.GeneDataCache;
import com.hartwig.hmftools.pave.ImpactClassifier;
import com.hartwig.hmftools.pave.PaveConfig;
import com.hartwig.hmftools.pave.VariantData;
import com.hartwig.hmftools.pave.VariantImpactBuilder;

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
    private final PaveConfig mConfig;
    private final ImpactClassifier mImpactClassifier;
    private final VariantImpactBuilder mImpactBuilder;
    private final GeneDataCache mGeneDataCache;

    private final List<String> mSampleIds;
    private final String mReferenceVariantsFile;
    private final DatabaseAccess mDbAccess;

    private final boolean mOnlyDriverGenes;
    private final boolean mOnlyCanonical;

    private BufferedWriter mCsvWriter;

    // counters
    private int mVariantCount;
    private int mTransVariantCount;
    private final PerformanceCounter mPerfCounter;

    private static final String SAMPLE_ID_FILE = "sample_id_file";
    private static final String REF_VARIANTS_FILE = "ref_variants_file";
    private static final String ONLY_DRIVER_GENES = "only_driver_genes";
    private static final String ONLY_CANONCIAL = "only_canonical";

    public ImpactComparisons(final CommandLine cmd)
    {
        mConfig = new PaveConfig(cmd);

        mSampleIds = Lists.newArrayList();

        mSampleIds.addAll(loadSampleIdsFile(cmd.getOptionValue(SAMPLE_ID_FILE)));

        mGeneDataCache = new GeneDataCache(
                cmd.getOptionValue(ENSEMBL_DATA_DIR), mConfig.RefGenVersion, cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION));

        mImpactBuilder = new VariantImpactBuilder(mGeneDataCache);

        RefGenomeInterface refGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));
        mImpactClassifier = new ImpactClassifier(refGenome);

        mDbAccess = createDatabaseAccess(cmd);
        mReferenceVariantsFile = cmd.getOptionValue(REF_VARIANTS_FILE);

        mOnlyCanonical = cmd.hasOption(ONLY_CANONCIAL);
        mOnlyDriverGenes = cmd.hasOption(ONLY_DRIVER_GENES);

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
        if(mSampleIds.isEmpty())
        {
            VI_LOGGER.error("missing sampleIds, exiting");
            System.exit(1);
        }

        if(mDbAccess == null && mReferenceVariantsFile == null)
        {
            VI_LOGGER.error("neither DB nor ref variants file configured, exiting");
            System.exit(1);
        }

        if(!mGeneDataCache.loadCache(mOnlyCanonical, mOnlyDriverGenes))
        {
            VI_LOGGER.error("Ensembl data cache loading failed, exiting");
            System.exit(1);
        }

        initialiseImpactWriter();

        if(mDbAccess != null)
        {
            for(int i = 0; i < mSampleIds.size(); ++i)
            {
                String sampleId = mSampleIds.get(i);

                mPerfCounter.start();
                loadSampleDatabaseRecords(sampleId);
                mPerfCounter.stop();

                if(i > 0 && (i % 100) == 0)
                {
                    VI_LOGGER.info("processing {} samples", i);
                }
            }
        }
        else
        {
            processRefVariantFile(mReferenceVariantsFile);
        }

        close();

        mPerfCounter.logStats();

        VI_LOGGER.info("impact comparison complete");
    }

    private void loadSampleDatabaseRecords(final String sampleId)
    {
        mVariantCount = 0;
        mTransVariantCount = 0;

        Result<Record> result = mDbAccess.context().select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(SOMATICVARIANT.GENE.isNotNull())
                .orderBy(SOMATICVARIANT.CHROMOSOME,SOMATICVARIANT.POSITION)
                .fetch();

        for (Record record : result)
        {
            final SomaticVariant variant = SomaticVariantDAO.buildFromRecord(record);
            processVariant(sampleId, RefVariantData.fromSomatic(variant));
        }

        VI_LOGGER.debug("sample({}) processed {} variants and transcripts({})",
                sampleId, mVariantCount, mTransVariantCount);
    }

    private void processRefVariantFile(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));
            String header = fileReader.readLine();

            final String fileDelim = "\t";
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, fileDelim);

            int sampleIndex = fieldsIndexMap.get("sampleId");

            int chrIndex = fieldsIndexMap.get("chromosome");
            int posIndex = fieldsIndexMap.get("position");
            int refIndex = fieldsIndexMap.get("ref");
            int altIndex = fieldsIndexMap.get("alt");
            int typeIndex = fieldsIndexMap.get("type");
            int geneIndex = fieldsIndexMap.get("gene");
            int canonicalEffectIndex = fieldsIndexMap.get("canonicalEffect");
            int canonicalCodingEffectIndex = fieldsIndexMap.get("canonicalCodingEffect");
            int worstEffectIndex = fieldsIndexMap.get("worstEffect");
            int worstCodingEffectIndex = fieldsIndexMap.get("worstCodingEffect");
            int worstEffectTranscriptIndex = fieldsIndexMap.get("worstEffectTranscript");
            int genesAffectedIndex = fieldsIndexMap.get("genesEffected");
            int canonicalHgvsCodingImpactIndex = fieldsIndexMap.get("canonicalHgvsCodingImpact");
            int canonicalHgvsProteinImpactIndex = fieldsIndexMap.get("canonicalHgvsProteinImpact");
            int microhomologyIndex = fieldsIndexMap.get("microhomology");
            int repeatSequenceIndex = fieldsIndexMap.get("repeatSequence");
            int phasedInframeIndelIndex = fieldsIndexMap.get("phasedInframeIndel");

            String currentSample = "";
            int samplesProcessed = 0;

            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(fileDelim, -1);

                String sampleId = items[sampleIndex];

                if(!currentSample.equals(sampleId))
                {
                    if(!mSampleIds.contains(sampleId))
                    {
                        if(samplesProcessed == mSampleIds.size())
                        {
                            VI_LOGGER.info("all samples processed");
                            return;
                        }
                        else
                        {
                            continue;
                        }
                    }

                    currentSample = sampleId;
                    ++samplesProcessed;

                    if(samplesProcessed > 0 && (samplesProcessed % 100) == 0)
                    {
                        VI_LOGGER.info("processing {} samples", samplesProcessed);
                    }
                }

                RefVariantData variant = new RefVariantData(
                        items[chrIndex], Integer.parseInt(items[posIndex]), items[refIndex], items[altIndex],
                        VariantType.valueOf(items[typeIndex]), items[geneIndex], items[canonicalEffectIndex],
                        items[canonicalCodingEffectIndex].isEmpty() ? NONE : CodingEffect.valueOf(items[canonicalCodingEffectIndex]),
                        items[worstEffectIndex], CodingEffect.valueOf(items[worstCodingEffectIndex]), items[worstEffectTranscriptIndex],
                        Integer.parseInt(items[genesAffectedIndex]), items[canonicalHgvsCodingImpactIndex], items[canonicalHgvsProteinImpactIndex],
                        items[microhomologyIndex], items[repeatSequenceIndex], Boolean.parseBoolean(items[phasedInframeIndelIndex]));

                processVariant(sampleId, variant);
            }
        }
        catch(IOException e)
        {
            VI_LOGGER.error("failed to read ref variant data file: {}", e.toString());
        }
    }

    private void processVariant(final String sampleId, final RefVariantData refVariant)
    {
        ++mVariantCount;

        // generate variant impact data and then write comparison results to CSV file
        VariantData variant = new VariantData(
                refVariant.Chromosome, refVariant.Position, refVariant.Ref, refVariant.Alt);

        variant.setVariantDetails(refVariant.PhasedInframeIndel, refVariant.Microhomology, refVariant.RepeatSequence);

        findVariantImpacts(variant, mImpactClassifier, mGeneDataCache);

        VariantImpact variantImpact = mImpactBuilder.createVariantImpact(variant);

        if(VI_LOGGER.isDebugEnabled())
        {
            if(hasCodingEffectDiff(variantImpact.CanonicalCodingEffect, refVariant.CanonicalCodingEffect))
            {
                VI_LOGGER.debug("sample({}) var({}) diff canonical coding: new({}) snpEff({})",
                        sampleId, variant.toString(), variantImpact.CanonicalCodingEffect, refVariant.CanonicalCodingEffect);
            }
        }

        writeVariantData(sampleId, variant, variantImpact, refVariant);
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
            VI_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return;
        }
    }

    private void writeVariantData(
            final String sampleId, final VariantData variant, final VariantImpact variantImpact, final RefVariantData refVariant)
    {
        try
        {
            mCsvWriter.write(String.format("%s,%s,%s,%s",
                    sampleId, variant.toCsv(), variantImpact.CanonicalGeneId, variantImpact.CanonicalGeneName));

            boolean isDriver = mGeneDataCache.getDriverPanelGenes().contains(variantImpact.CanonicalGeneName);

            if(!isDriver && !variantImpact.CanonicalGeneName.equals(refVariant.Gene))
                isDriver = mGeneDataCache.getDriverPanelGenes().contains(refVariant.Gene);

            mCsvWriter.write(String.format(",%s,%s,%s,%s,%s,%s,%d",
                    isDriver, variantImpact.CanonicalEffect, variantImpact.CanonicalCodingEffect,
                    variantImpact.WorstEffect, variantImpact.WorstCodingEffect, variantImpact.WorstTranscript, variantImpact.GenesAffected));

            mCsvWriter.write(String.format(",%s,%s,%s,%s,%s,%s,%d",
                    refVariant.Gene, refVariant.CanonicalEffect, refVariant.CanonicalCodingEffect, refVariant.WorstEffect,
                    refVariant.WorstCodingEffect, refVariant.WorstEffectTranscript, refVariant.GenesAffected));

            mCsvWriter.newLine();
        }
        catch(IOException e)
        {
            VI_LOGGER.error("failed to write variant CSV file: {}", e.toString());
            return;
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = PaveConfig.createOptions();

        addDatabaseCmdLineArgs(options);
        options.addOption(SAMPLE_ID_FILE, true, "Sample ID file");
        options.addOption(REF_VARIANTS_FILE, true, "File with variants to test against");
        options.addOption(ONLY_DRIVER_GENES, false, "Only compare variants in driver genes");
        options.addOption(ONLY_CANONCIAL, false, "Only compare variants by canonical transcripts");

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
