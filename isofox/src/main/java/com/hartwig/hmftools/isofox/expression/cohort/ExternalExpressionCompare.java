package com.hartwig.hmftools.isofox.expression.cohort;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.GENE_DISTRIBUTION;
import static com.hartwig.hmftools.isofox.cohort.CohortAnalysisType.TRANSCRIPT_DISTRIBUTION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig.EXT_SOURCE_RSEM;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig.EXT_SOURCE_SALMON;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionCohortConfig.SOURCE_ISOFOX;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionData.fromIsofoxGene;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionData.fromIsofoxTranscript;
import static com.hartwig.hmftools.isofox.expression.cohort.ExpressionData.fromSalmon;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_ID;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_GENE_NAME;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.FLD_TRANS_ID;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_FITTED_FRAGMENTS;
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
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
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
        
        mGeneTransCache = new EnsemblDataCache(config.EnsemblDataCache, RefGenomeVersion.HG37);
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
                final String geneName = mGeneTransCache.getGeneDataByName(geneId).GeneName;

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
                final String geneName = mGeneTransCache.getGeneDataByName(geneId).GeneName;

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

        if(!formSampleFilenames(mConfig, mTransScope ? TRANSCRIPT_DISTRIBUTION : GENE_DISTRIBUTION, filenames))
            return;

        initialiseWriter();

        // load each sample's expression data, normalise and compare and the write results
        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path transcriptsFile = filenames.get(i);

            final Map<String,ExpressionData> isofoData = loadIsofoxFile(transcriptsFile);
            final Map<String,ExpressionData> externalData = loadFile(sampleId, mConfig.Expression.ExternalSource);

            ISF_LOGGER.debug("{}: sample({}) loaded transcript data", i, sampleId);
        }

        ISF_LOGGER.info("loaded {} samples transcript files", mConfig.SampleData.SampleIds.size());


        closeBufferedWriter(mWriter);
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mTransScope ? 
                    mConfig.formCohortFilename("trans_exp_compare.csv") : mConfig.formCohortFilename("gene_exp_compare.csv");
            
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("GeneId,GeneName,TransName");

            if(mTransScope)
                mWriter.write(",TransName");

            mWriter.write(String.format(",ISOFOX_TPM,%s_TPM,", mConfig.Expression.ExternalSource));

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write transcript data file: {}", e.toString());
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
            int transIdIndex = fieldsMap.get(FLD_TRANS_ID);
            int readCountIndex = fieldsMap.get(FLD_FITTED_FRAGMENTS);
            int tpmIndex = fieldsMap.get(FLD_TPM);

            for(final String data : lines)
            {
                ExpressionData expData = mTransScope ?
                        fromIsofoxTranscript(data, geneIdIndex, geneNameIndex, transIdIndex, readCountIndex, tpmIndex) :
                        fromIsofoxGene(data, geneIdIndex, geneNameIndex, readCountIndex, tpmIndex);

                if(expData == null)
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
        final String filename = source.equals(EXT_SOURCE_SALMON) ? String.format("%s.salmon.tsv", sampleId) : null;

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
                }
                if(source.equals(EXT_SOURCE_RSEM))
                {
                    // expData = fromSalmon(data, mGeneTransMap);
                }

                if(expData == null)
                    continue;

                if(mTransScope)
                {
                    expressionDataMap.put(expData.TransName, expData);
                }
                else
                {
                    ExpressionData geneExpData = expressionDataMap.get(expData.GeneId);

                    if(geneExpData == null)
                    {
                        geneExpData = new ExpressionData(source, expData.GeneId, expData.GeneName, "", expData.readCount(), expData.tpm());
                        expressionDataMap.put(expData.GeneId, geneExpData);
                    }
                    else
                    {
                        geneExpData.addCounts(expData.tpm(), expData.readCount());
                    }
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load external expression data file({}): {}", filename.toString(), e.toString());
        }

        return expressionDataMap;
    }
}
