package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanelConfig.DRIVER_GENE_PANEL_OPTION;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.variant.CodingEffect.NONE;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.addDatabaseCmdLineArgs;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.tables.Somaticvariant.SOMATICVARIANT;
import static com.hartwig.hmftools.pave.PaveApplication.findVariantImpacts;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;
import static com.hartwig.hmftools.pave.VariantData.NO_LOCAL_PHASE_SET;
import static com.hartwig.hmftools.pave.compare.ComparisonUtils.hasCodingEffectDiff;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

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
import com.hartwig.hmftools.pave.PhasedVariantClassifier;
import com.hartwig.hmftools.pave.PhasedVariants;
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
    private final ComparisonConfig mConfig;
    private final ImpactClassifier mImpactClassifier;
    private final VariantImpactBuilder mImpactBuilder;
    private final GeneDataCache mGeneDataCache;

    private final DatabaseAccess mDbAccess;

    private final ComparisonWriter mWriter;

    // counters
    private int mTotalComparisons;
    private int mMatchedCount;
    private final PerformanceCounter mPerfCounter;

    public ImpactComparisons(final CommandLine cmd)
    {
        mConfig = new ComparisonConfig(cmd);

        mGeneDataCache = new GeneDataCache(
                cmd.getOptionValue(ENSEMBL_DATA_DIR), mConfig.RefGenVersion, cmd.getOptionValue(DRIVER_GENE_PANEL_OPTION),
                true, false);

        mImpactBuilder = new VariantImpactBuilder(mGeneDataCache);

        RefGenomeInterface refGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));
        mImpactClassifier = new ImpactClassifier(refGenome);

        mDbAccess = createDatabaseAccess(cmd);

        mWriter = new ComparisonWriter(mGeneDataCache, mConfig);

        mPerfCounter = new PerformanceCounter("ClassifyVariant");
        mTotalComparisons = 0;
        mMatchedCount = 0;
    }

    public void run()
    {
        if(mConfig.SampleIds.isEmpty() && mConfig.ReferenceVariantsFile == null)
        {
            PV_LOGGER.error("missing sampleIds, exiting");
            System.exit(1);
        }

        if(mDbAccess == null && mConfig.ReferenceVariantsFile == null)
        {
            PV_LOGGER.error("neither DB nor ref variants file configured, exiting");
            System.exit(1);
        }

        if(!mGeneDataCache.loadCache(mConfig.OnlyCanonical, mConfig.OnlyDriverGenes))
        {
            PV_LOGGER.error("Ensembl data cache loading failed, exiting");
            System.exit(1);
        }

        if(mDbAccess != null)
        {
            for(int i = 0; i < mConfig.SampleIds.size(); ++i)
            {
                String sampleId = mConfig.SampleIds.get(i);

                mPerfCounter.start();
                loadSampleDatabaseRecords(sampleId);
                mPerfCounter.stop();

                if(i > 0 && (i % 100) == 0)
                {
                    PV_LOGGER.info("processing {} samples", i);
                }
            }
        }
        else
        {
            processRefVariantFile(mConfig.ReferenceVariantsFile);
        }

        mWriter.close();

        PV_LOGGER.info("total comparisons({}) matched({}) diffs({})",
                mTotalComparisons, mMatchedCount, mTotalComparisons - mMatchedCount);

        mPerfCounter.logStats();

        PV_LOGGER.info("impact comparison complete");
    }

    private void processVariant(final String sampleId, final RefVariantData refVariant)
    {
        // generate variant impact data and then write comparison results to CSV file
        VariantData variant = new VariantData(
                refVariant.Chromosome, refVariant.Position, refVariant.Ref, refVariant.Alt);

        variant.setVariantDetails(refVariant.LocalPhaseSet, refVariant.Microhomology, refVariant.RepeatSequence);
        variant.setSampleId(sampleId);
        variant.setRefData(refVariant);

        findVariantImpacts(variant, mImpactClassifier, mGeneDataCache);

        processPhasedVariants(variant.localPhaseSet());

        if(!variant.hasLocalPhaseSet())
            processVariant(sampleId, variant, refVariant);
    }

    private void processPhasedVariants(int currentLocalPhaseSet)
    {
        List<VariantData> variants = mImpactClassifier.processPhasedVariants(currentLocalPhaseSet);

        if(variants != null)
            variants.forEach(x -> processVariant(x.sampleId(), x, x.refData()));
    }

    private void processVariant(final String sampleId, final VariantData variant, final RefVariantData refVariant)
    {
        ++mTotalComparisons;

        VariantImpact variantImpact = mImpactBuilder.createVariantImpact(variant);

        if(!hasCodingEffectDiff(variantImpact.CanonicalCodingEffect, refVariant.CanonicalCodingEffect))
        {
            ++mMatchedCount;
            return;
        }

        logComparison(sampleId, refVariant, variant, variantImpact);

        mWriter.writeVariantData(sampleId, variant, variantImpact, refVariant);
    }

    private void logComparison(
            final String sampleId, final RefVariantData refVariant, final VariantData variant, final VariantImpact variantImpact)
    {
        if(!PV_LOGGER.isDebugEnabled())
            return;

        if(variant.getImpacts().isEmpty())
        {
            PV_LOGGER.debug("sample({}) var({}) no gene impacts found vs snpEff({} : {})",
                    sampleId, variant.toString(), refVariant.Gene, refVariant.CanonicalCodingEffect);
            return;
        }

        if(!variant.getImpacts().containsKey(refVariant.Gene))
        {
            PV_LOGGER.debug("sample({}) var({}) diff gene and canonical coding: pave({} : {}) snpEff({} : {})",
                    sampleId, variant.toString(), variantImpact.CanonicalGeneName, variantImpact.CanonicalCodingEffect,
                    refVariant.Gene, refVariant.CanonicalCodingEffect);
            return;
        }

        PV_LOGGER.debug("sample({}) var({}) diff canonical coding: pave({} : {}) snpEff({} : {})",
                sampleId, variant.toString(), variantImpact.CanonicalGeneName, variantImpact.CanonicalCodingEffect,
                refVariant.Gene, refVariant.CanonicalCodingEffect);
    }

    private void loadSampleDatabaseRecords(final String sampleId)
    {
        mTotalComparisons = 0;
        mMatchedCount = 0;

        Result<Record> result = mDbAccess.context().select()
                .from(SOMATICVARIANT)
                .where(SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(SOMATICVARIANT.GENE.isNotNull())
                .orderBy(SOMATICVARIANT.CHROMOSOME,SOMATICVARIANT.POSITION)
                .fetch();

        for(Record record : result)
        {
            final SomaticVariant variant = SomaticVariantDAO.buildFromRecord(record);
            processVariant(sampleId, RefVariantData.fromSomatic(variant));
        }

        processPhasedVariants(NO_LOCAL_PHASE_SET);
        mImpactClassifier.phasedVariants().clear();

        PV_LOGGER.debug("sample({}) processed {} variants and transcripts({})",
                sampleId, mTotalComparisons, mMatchedCount);
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
            int worstCodingEffectIndex = fieldsIndexMap.get("worstCodingEffect");
            int genesAffectedIndex = fieldsIndexMap.get("genesEffected");
            int canonicalHgvsCodingImpactIndex = fieldsIndexMap.get("canonicalHgvsCodingImpact");
            int canonicalHgvsProteinImpactIndex = fieldsIndexMap.get("canonicalHgvsProteinImpact");
            int microhomologyIndex = fieldsIndexMap.get("microhomology");
            int repeatSequenceIndex = fieldsIndexMap.get("repeatSequence");
            int phasedInframeIndelIndex = fieldsIndexMap.get("phasedInframeIndel");
            int localPhaseSetIndex = fieldsIndexMap.get("localPhaseSet");
            int reportedIndex = fieldsIndexMap.get("reported");

            String currentSample = "";
            int samplesProcessed = 0;

            String line = "";

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(fileDelim, -1);

                String sampleId = items[sampleIndex];

                if(!currentSample.equals(sampleId))
                {
                    if(!mConfig.SampleIds.isEmpty() && !mConfig.SampleIds.contains(sampleId))
                    {
                        if(samplesProcessed == mConfig.SampleIds.size())
                        {
                            PV_LOGGER.info("all {} samples processed", samplesProcessed);
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
                        PV_LOGGER.info("processed {} samples", samplesProcessed);
                    }
                }

                int localPhaseSet = !items[localPhaseSetIndex].equals("NULL") ? Integer.parseInt(items[localPhaseSetIndex]) : NO_LOCAL_PHASE_SET;

                RefVariantData variant = new RefVariantData(
                        items[chrIndex], Integer.parseInt(items[posIndex]), items[refIndex], items[altIndex],
                        VariantType.valueOf(items[typeIndex]), items[geneIndex], items[canonicalEffectIndex],
                        items[canonicalCodingEffectIndex].isEmpty() ? NONE : CodingEffect.valueOf(items[canonicalCodingEffectIndex]),
                        CodingEffect.valueOf(items[worstCodingEffectIndex]), Integer.parseInt(items[genesAffectedIndex]),
                        items[canonicalHgvsCodingImpactIndex], items[canonicalHgvsProteinImpactIndex],
                        items[microhomologyIndex], items[repeatSequenceIndex], Boolean.parseBoolean(items[phasedInframeIndelIndex]),
                        localPhaseSet, Boolean.parseBoolean(items[reportedIndex]));

                processVariant(sampleId, variant);
            }

            processPhasedVariants(NO_LOCAL_PHASE_SET);
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to read ref variant data file: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = ComparisonConfig.createOptions();

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
