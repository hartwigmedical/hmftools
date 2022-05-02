package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_TRANS_NAME;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.GENE_EXPRESSION_MATRIX;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.TRANSCRIPT_EXPRESSION_MATRIX;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig.EXT_SOURCE_RSEM;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig.EXT_SOURCE_SALMON;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionData.fromIsofoxGene;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionData.fromIsofoxTranscript;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionData.fromRsemGene;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionData.fromRsemTranscript;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionData.fromSalmon;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionData.getExternalSourceFilename;
import static com.hartwig.hmftools.isofox.results.GeneResult.FLD_SPLICED_FRAGS;
import static com.hartwig.hmftools.isofox.results.GeneResult.FLD_UNSPLICED_FRAGS;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_EFFECTIVE_LENGTH;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_FITTED_FRAGMENTS;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_RAW_FRAGMENTS;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class ExternalExpressionCompare
{
    private final CohortConfig mConfig;

    private final EnsemblDataCache mGeneTransCache;
    private final Map<String,String[]> mGeneTransMap; // geneId and geneName keyed by transName
    private final boolean mTransScope;

    private BufferedWriter mWriter;

    public ExternalExpressionCompare(final CohortConfig config)
    {
        mConfig = config;
        mTransScope = config.Expression.TranscriptScope;
        
        mGeneTransMap = Maps.newHashMap();
        
        mGeneTransCache = new EnsemblDataCache(config.EnsemblDataCache, RefGenomeVersion.V37);
        mGeneTransCache.setRequiredData(false, false, false, false);
        buildTransGeneMap();
    }

    private void buildTransGeneMap()
    {
        if(mConfig.RestrictedGeneIds.isEmpty())
        {
            mGeneTransCache.load(false);

            for(Map.Entry<String,List<TranscriptData>> entry : mGeneTransCache.getTranscriptDataMap().entrySet())
            {
                final String geneId = entry.getKey();
                final String geneName = mGeneTransCache.getGeneDataById(geneId).GeneName;

                for(TranscriptData transData : entry.getValue())
                {
                    mGeneTransMap.put(transData.TransName, new String[] { geneId, geneName});
                }
            }
        }
        else
        {
            mGeneTransCache.load(true);
            mGeneTransCache.loadTranscriptData(mConfig.RestrictedGeneIds);

            for(final String geneId : mConfig.RestrictedGeneIds)
            {
                final String geneName = mGeneTransCache.getGeneDataById(geneId).GeneName;

                for(TranscriptData transData : mGeneTransCache.getTranscripts(geneId))
                {
                    mGeneTransMap.put(transData.TransName, new String[] { geneId, geneName});
                }
            }
        }
    }

    public void processSampleTranscriptFiles()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, mTransScope ? TRANSCRIPT_EXPRESSION_MATRIX : GENE_EXPRESSION_MATRIX, filenames))
            return;

        initialiseWriter();

        // load each sample's expression data, normalise and compare and the write results
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path transcriptsFile = filenames.get(i);

            final Map<String,ExpressionData> isofoxExpMap = loadIsofoxFile(transcriptsFile);
            final Map<String,ExpressionData> externalExpMap = loadFile(sampleId, mConfig.Expression.ExternalSource);

            ISF_LOGGER.debug("{}: sample({}) loaded expression data", i, sampleId);

            processExpressionData(sampleId, isofoxExpMap, externalExpMap);
        }

        ISF_LOGGER.info("processed {} samples transcript files", mConfig.SampleData.SampleIds.size());

        closeBufferedWriter(mWriter);
    }

    private void processExpressionData(final String sampleId, final Map<String,ExpressionData> isofoxExpMap, final Map<String,ExpressionData> externalExpMap)
    {
        for(Map.Entry<String,ExpressionData> entry : isofoxExpMap.entrySet())
        {
            final String id = entry.getKey();
            final ExpressionData isofoxExpData = entry.getValue();
            final ExpressionData externalExpData = externalExpMap.get(id);

            // apply any filters..

            writeComparisonData(sampleId, isofoxExpData, externalExpData);
        }
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mTransScope ?
                    mConfig.formCohortFilename("trans_exp_compare.csv") : mConfig.formCohortFilename("gene_exp_compare.csv");

            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("SampleId,GeneId,GeneName");

            final String sourceName = mConfig.Expression.ExternalSource;

            if(mTransScope)
            {
                mWriter.write(String.format(",TransName,ISOFOX_TPM,%s_TPM,RawFrags,FittedFrags,%s_Reads,EffectiveLength,LowMapQualFrags",
                        sourceName, sourceName));
            }
            else
            {
                mWriter.write(String.format(",ISOFOX_TPM,%s_TPM,SplicedFragments,UnsplicedFragments,LowMapQualFrags", sourceName));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write expression comparison data file: {}", e.toString());
        }
    }

    private void writeComparisonData(final String sampleId, final ExpressionData isofoxExpData, final ExpressionData externalExpData)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s", sampleId, isofoxExpData.GeneId, isofoxExpData.GeneName));

            if(mTransScope)
            {
                mWriter.write(String.format(",%s,%.3g,%.3g,%.1f,%.1f,%d,%d,%.1f",
                        isofoxExpData.TransName, isofoxExpData.tpm(), externalExpData != null ? externalExpData.tpm() : -1,
                        isofoxExpData.rawFragment(), isofoxExpData.fittedFragments(),
                        externalExpData != null ? externalExpData.readCount() : -1, isofoxExpData.EffectiveLength, isofoxExpData.LowMapQualFrags));
            }
            else
            {
                mWriter.write(String.format(",%.3g,%.3g,%d,%d,%.1f",
                        isofoxExpData.tpm(), externalExpData != null ? externalExpData.tpm() : -1,
                        isofoxExpData.SplicedFragments, isofoxExpData.UnsplicedFragments, isofoxExpData.LowMapQualFrags));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write expression comparison data file: {}", e.toString());
        }
    }

    private Map<String,ExpressionData> loadIsofoxFile(final Path filename)
    {
        final Map<String,ExpressionData> expressionDataMap = Maps.newHashMap();

        try
        {
            final List<String> lines = Files.readAllLines(filename);

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            int geneIdIndex = fieldsMap.get(FLD_GENE_ID);
            int geneNameIndex = fieldsMap.get(FLD_GENE_NAME);
            int tpmIndex = fieldsMap.get(FLD_TPM);

            Integer fittedFragsIndex = fieldsMap.get(FLD_FITTED_FRAGMENTS);
            Integer rawFragsIndex = fieldsMap.get(FLD_RAW_FRAGMENTS);
            Integer transNameIndex = fieldsMap.get(FLD_TRANS_NAME);
            Integer effectiveLengthIndex = fieldsMap.get(FLD_EFFECTIVE_LENGTH);
            Integer splicedIndex = fieldsMap.get(FLD_SPLICED_FRAGS);
            Integer unsplicedIndex = fieldsMap.get(FLD_UNSPLICED_FRAGS);
            Integer lowQualIndex = fieldsMap.get("LowMapQualFrags");

            for(final String data : lines)
            {
                ExpressionData expData = mTransScope ?
                        fromIsofoxTranscript(
                                data, geneIdIndex, geneNameIndex, transNameIndex, fittedFragsIndex, rawFragsIndex,
                                tpmIndex, effectiveLengthIndex, lowQualIndex) :
                        fromIsofoxGene(data, geneIdIndex, geneNameIndex, tpmIndex, splicedIndex, unsplicedIndex, lowQualIndex);

                if(expData == null)
                    continue;

                if(!mConfig.RestrictedGeneIds.isEmpty() && !mConfig.RestrictedGeneIds.contains(expData.GeneId))
                    continue;

                if(mTransScope)
                    expressionDataMap.put(expData.TransName, expData);
                else
                    expressionDataMap.put(expData.GeneId, expData);
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load Isofox expression data file({}): {}", filename.toString(), e.toString());
        }

        return expressionDataMap;
    }

    private Map<String,ExpressionData> loadFile(final String sampleId, final String source)
    {
        final Map<String,ExpressionData> expressionDataMap = Maps.newHashMap();

        final String filename = String.format("%s/%s", mConfig.RootDataDir, getExternalSourceFilename(source, sampleId, mTransScope));

        if(filename == null)
            return expressionDataMap;

        try
        {

            final List<String> lines = Files.readAllLines(Paths.get(filename));
            lines.remove(0);

            for(final String data : lines)
            {
                ExpressionData expData = null;

                if(source.equals(EXT_SOURCE_SALMON))
                {
                    expData = fromSalmon(data, mGeneTransMap);

                    if(expData == null)
                        continue;

                    if(mTransScope)
                    {
                        expressionDataMap.put(expData.TransName, expData);
                    }
                    else
                    {
                        // build up gene data from transcripts
                        ExpressionData geneExpData = expressionDataMap.get(expData.GeneId);

                        if(geneExpData == null)
                        {
                            geneExpData = new ExpressionData(source, expData.GeneId, expData.GeneName, "",
                                    0, 0, expData.readCount(), expData.tpm(),
                                    0, 0, 0, 0);
                            expressionDataMap.put(expData.GeneId, geneExpData);
                        }
                        else
                        {
                            geneExpData.addCounts(expData.tpm(), 0, 0, expData.readCount());
                        }
                    }
                }
                else if(source.equals(EXT_SOURCE_RSEM))
                {
                    if(mTransScope)
                        expData = fromRsemTranscript(data, mGeneTransMap);
                    else
                        expData = fromRsemGene(data, mGeneTransCache);

                    if(expData == null)
                        continue;

                    if(mTransScope)
                        expressionDataMap.put(expData.TransName, expData);
                    else
                        expressionDataMap.put(expData.GeneId, expData);
                }

            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load external expression data file({}): {}", filename, e.toString());
        }

        return expressionDataMap;
    }
}
