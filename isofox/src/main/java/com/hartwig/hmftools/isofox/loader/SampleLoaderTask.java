package com.hartwig.hmftools.isofox.loader;

import static com.hartwig.hmftools.common.utils.config.ConfigUtils.convertWildcardSamplePath;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.expression.cohort.CohortGenePercentiles.CANCER_TYPE_OTHER;
import static com.hartwig.hmftools.isofox.loader.DataLoadType.NOVEL_JUNCTION;
import static com.hartwig.hmftools.patientdb.dao.DatabaseAccess.createDatabaseAccess;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.Callable;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.GeneExpressionFile;
import com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.RnaFusionFile;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaStatisticFile;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
import com.hartwig.hmftools.patientdb.dao.IsofoxDAO;

import org.jooq.Record;
import org.jooq.Result;

public class SampleLoaderTask implements Callable<Void>
{
    private final int mTaskId;
    private final DataLoaderConfig mConfig;

    private final List<String> mSampleIds;
    private final DatabaseAccess mDbAccess;
    private final IsofoxDAO mRnaDAO;

    public SampleLoaderTask(int taskId, final DataLoaderConfig config)
    {
        mTaskId = taskId;
        mConfig = config;

        mSampleIds = Lists.newArrayList();

        mDbAccess = createDatabaseAccess(mConfig.ConfigItems);

        if(mDbAccess == null)
        {
            ISF_LOGGER.error("invalid DB connection");
            System.exit(1);
            mRnaDAO = null;
        }
        else
        {
            mRnaDAO = new IsofoxDAO(mDbAccess.context());
        }
    }

    public final List<String> getSampleIds() { return mSampleIds; }

    @Override
    public Void call()
    {
        if(mRnaDAO == null)
            return null;

        ISF_LOGGER.info("{}: loading data for {} samples", mTaskId, mSampleIds.size());

        for(int i = 0; i < mSampleIds.size(); ++i)
        {
            String sampleId = mSampleIds.get(i);

            loadSampleData(sampleId);

            if(i > 0 && (i % 100) == 0)
            {
                ISF_LOGGER.info("{}: processed {} samples", mTaskId, i);
            }
        }

        if(mConfig.Threads > 1)
        {
            ISF_LOGGER.info("{}: tasks complete for {} samples", mTaskId, mSampleIds.size());
        }

        return null;
    }

    private void loadSampleData(final String sampleId)
    {
        final String cancerType = getSampleCancerType(sampleId);

        if(cancerType == null)
        {
            ISF_LOGGER.warn("sample({}) data loading skipped", sampleId);
            return;
        }

        ISF_LOGGER.debug("sample({}) cancerType({}) loading Isofox RNA data", sampleId, cancerType);

        loadStatistics(sampleId);
        loadGeneExpression(sampleId, cancerType);
        loadNovelJunctions(sampleId);
        loadFusions(sampleId);
    }

    private String getSampleCancerType(final String sampleId)
    {
        if(mConfig.SampleCancerTypes.containsKey(sampleId))
            return mConfig.SampleCancerTypes.get(sampleId);

        try
        {
            final String queryStr = String.format("select primaryTumorLocation from clinical where sampleId = '%s';", sampleId);

            final Result<Record> result = mDbAccess.context().fetch(queryStr);

            String primaryTumorLocation = CANCER_TYPE_OTHER;

            for(Record record : result)
            {
                Object ptLocation = record.getValue("primaryTumorLocation");
                primaryTumorLocation = ptLocation != null ? ptLocation.toString() : CANCER_TYPE_OTHER;
                break;
            }

            if(mConfig.PrimaryCancerTypes.contains(primaryTumorLocation))
                return primaryTumorLocation;
        }
        catch(Exception e)
        {
            ISF_LOGGER.warn("failed to retrieve clinical data for sample({}): {}", sampleId, e.toString());
            return null;
        }

        return CANCER_TYPE_OTHER;
    }

    private void loadStatistics(final String sampleId)
    {
        if(!mConfig.loadDataType(DataLoadType.STATISTICS))
            return;

        String sampleDataDir = convertWildcardSamplePath(mConfig.StatisticsDataDir, sampleId);

        String filename = RnaStatisticFile.generateFilename(sampleDataDir, sampleId);

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            if(lines.size() != 2)
            {
                ISF_LOGGER.error("sample({}) invalid summary file", sampleId);
                return;
            }

            RnaStatistics statistics = RnaStatisticFile.fromLines(lines);

            ISF_LOGGER.debug("sample({}) writing summary statistics to DB", sampleId);
            mRnaDAO.writeRnaStatistics(sampleId, statistics);

        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load summary statistics data file({}): {}", filename, e.toString());
        }
    }

    private void loadGeneExpression(final String sampleId, final String cancerType)
    {
        if(!mConfig.loadDataType(DataLoadType.GENE_EXPRESSION))
            return;

        String sampleDataDir = convertWildcardSamplePath(mConfig.GeneDataDir, sampleId);

        String filename = GeneExpressionFile.generateFilename(sampleDataDir, sampleId);

        List<GeneExpression> geneExpressions = GeneExpressionFile.read(filename);

        ISF_LOGGER.debug("sample({}) writing {} gene expression records to DB", sampleId, geneExpressions.size());
        mRnaDAO.writeGeneExpressions(sampleId, geneExpressions);
    }

    private void loadNovelJunctions(final String sampleId)
    {
        if(!mConfig.loadDataType(NOVEL_JUNCTION))
            return;

        String sampleDataDir = convertWildcardSamplePath(mConfig.AltSjDataDir, sampleId);

        String filename = NovelSpliceJunctionFile.generateFilename(sampleDataDir, sampleId);
        List<NovelSpliceJunction> novelJunctions = NovelSpliceJunctionFile.read(filename);

        // no restricted gene lists any more
        /*
        if(!mConfig.processGeneId(geneId))
            continue;
        */

        ISF_LOGGER.debug("sample({}) writing {} novel SJs records to DB", sampleId, novelJunctions.size());
        mRnaDAO.writeNovelSpliceJunctions(sampleId, novelJunctions);
    }

    private void loadFusions(final String sampleId)
    {
        if(!mConfig.loadDataType(DataLoadType.FUSION))
            return;

        String sampleDataDir = convertWildcardSamplePath(mConfig.FusionDataDir, sampleId);

        String filename = RnaFusionFile.generateFilename(sampleDataDir, sampleId);

        List<RnaFusion> fusions = RnaFusionFile.read(filename);

        ISF_LOGGER.debug("sample({}) writing {} fusion records to DB", sampleId, fusions.size());
        mRnaDAO.writeRnaFusions(sampleId, fusions);
    }

}
