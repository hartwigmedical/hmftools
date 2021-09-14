package com.hartwig.hmftools.isofox.loader;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.ALT_SJ_FILE_ID;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_FRAG_COUNT;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_TYPE;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.formKey;
import static com.hartwig.hmftools.common.rna.GeneExpressionFile.GENE_DATA_FILE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_ID;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_GENE_NAME;
import static com.hartwig.hmftools.common.rna.RnaCommon.ISF_FILE_ID;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.expression.cohort.SampleGenePercentiles.CANCER_TYPE_OTHER;
import static com.hartwig.hmftools.isofox.expression.cohort.SampleGenePercentiles.PAN_CANCER;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_CHR;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_COHORT_COUNT;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_COVERAGE;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_DISCORD_FRAGS;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_JUNC_TYPE;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_MAX_ANCHOR;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_ORIENT;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_POS;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_REALIGN_FLAGS;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_SPLIT_FRAGS;
import static com.hartwig.hmftools.isofox.fusion.FusionData.FLD_SV_TYPE;
import static com.hartwig.hmftools.isofox.fusion.FusionData.formStreamField;
import static com.hartwig.hmftools.isofox.fusion.FusionWriter.PASS_FUSION_FILE_ID;
import static com.hartwig.hmftools.isofox.loader.DataLoadType.NOVEL_JUNCTION;
import static com.hartwig.hmftools.isofox.results.GeneResult.FLD_SPLICED_FRAGS;
import static com.hartwig.hmftools.isofox.results.GeneResult.FLD_UNSPLICED_FRAGS;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.SUMMARY_FILE;
import static com.hartwig.hmftools.isofox.results.TranscriptResult.FLD_TPM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.ImmutableGeneExpression;
import com.hartwig.hmftools.common.rna.ImmutableNovelSpliceJunction;
import com.hartwig.hmftools.common.rna.ImmutableRnaFusion;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.common.rna.RnaStatistics;
import com.hartwig.hmftools.isofox.expression.cohort.SampleGenePercentiles;
import com.hartwig.hmftools.patientdb.dao.IsofoxDAO;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;
import org.jooq.Record;
import org.jooq.Result;

public class IsofoxDataLoader
{
    private final DataLoaderConfig mConfig;
    private final SampleGenePercentiles mGeneDistribution;
    private final Map<String,Integer> mAltSjCohortFrequency;

    public IsofoxDataLoader(final CommandLine cmd)
    {
        mConfig = new DataLoaderConfig(cmd);
        mGeneDistribution = new SampleGenePercentiles(mConfig.GeneDistributionFile);

        mAltSjCohortFrequency = Maps.newHashMap();
        loadAltSjCohortFile();
    }

    public boolean load()
    {
        if(mConfig.SampleIds.isEmpty())
        {
            ISF_LOGGER.error("no sample IDs configured");
            return false;
        }

        if(mConfig.DbAccess == null)
        {
            ISF_LOGGER.error("invalid DB connection");
            return false;
        }

        final IsofoxDAO rnaDAO = new IsofoxDAO(mConfig.DbAccess.context());

        for(String sampleId : mConfig.SampleIds)
        {
            final String cancerType = getSampleCancerType(sampleId);

            ISF_LOGGER.info("sample({}) cancerType({}) loading Isofox RNA data", sampleId, cancerType);

            loadStatistics(sampleId, rnaDAO);
            loadGeneExpression(sampleId, cancerType, rnaDAO);
            loadNovelJunctions(sampleId, rnaDAO);
            loadFusions(sampleId, rnaDAO);
        }

        return true;
    }

    private String getSampleCancerType(final String sampleId)
    {
        if(mConfig.SampleCancerTypes.containsKey(sampleId))
            return mConfig.SampleCancerTypes.get(sampleId);

        try
        {
            final String queryStr = String.format("select primaryTumorLocation from clinical where sampleId = '%s';", sampleId);

            final Result<Record> result = mConfig.DbAccess.context().fetch(queryStr);

            String primaryTumorLocation = CANCER_TYPE_OTHER;

            for(Record record : result)
            {
                primaryTumorLocation = record.getValue("primaryTumorLocation").toString();
                break;
            }

            if(mConfig.PrimaryCancerTypes.contains(primaryTumorLocation))
                return primaryTumorLocation;
        }
        catch(Exception e)
        {
            ISF_LOGGER.error("failed to retrieve clinical data for sample({}): {}", e.toString());
        }

        return CANCER_TYPE_OTHER;
    }

    private void loadStatistics(final String sampleId, final IsofoxDAO rnaDAO)
    {
        if(!mConfig.loadDataType(DataLoadType.STATISTICS))
            return;

        final String sampleDataDir = mConfig.StatisticsDataDir.contains("*") ?
                mConfig.StatisticsDataDir.replaceAll("\\*", sampleId) : mConfig.StatisticsDataDir;

        final String filename = sampleDataDir + sampleId + ISF_FILE_ID + SUMMARY_FILE;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            if(lines.size() != 2)
            {
                ISF_LOGGER.error("sample({}) invalid summary file", sampleId);
                return;
            }

            final RnaStatistics statistics = RnaStatistics.fromCsv(lines.get(1));

            ISF_LOGGER.debug("sample({}) writing summary statistics to DB", sampleId);
            rnaDAO.writeRnaStatistics(sampleId, statistics);

        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load summary statistics data file({}): {}", filename, e.toString());
            return;
        }
    }

    private void loadGeneExpression(final String sampleId, final String cancerType, final IsofoxDAO rnaDAO)
    {
        if(!mConfig.loadDataType(DataLoadType.GENE_EXPRESSION))
            return;

        final List<GeneExpression> geneExpressions = Lists.newArrayList();

        final String sampleDataDir = mConfig.GeneDataDir.contains("*") ?
                mConfig.GeneDataDir.replaceAll("\\*", sampleId) : mConfig.GeneDataDir;

        final String filename = sampleDataDir + sampleId + ISF_FILE_ID + GENE_DATA_FILE_ID;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                final String geneId = items[fieldsIndexMap.get(FLD_GENE_ID)];

                if(!mConfig.processGeneId(geneId))
                    continue;

                double tpm = Double.parseDouble(items[fieldsIndexMap.get(FLD_TPM)]);

                double medianCancer = mGeneDistribution.getTpmMedian(geneId, cancerType);
                double percentileCancer = mGeneDistribution.getTpmPercentile(geneId, cancerType, tpm);
                double medianCohort = mGeneDistribution.getTpmMedian(geneId, PAN_CANCER);
                double percentileCohort = mGeneDistribution.getTpmPercentile(geneId, PAN_CANCER, tpm);

                geneExpressions.add(ImmutableGeneExpression.builder()
                        .geneName(items[fieldsIndexMap.get(FLD_GENE_NAME)])
                        .tpm(tpm)
                        .splicedFragments(Integer.parseInt(items[fieldsIndexMap.get(FLD_SPLICED_FRAGS)]))
                        .unsplicedFragments(Integer.parseInt(items[fieldsIndexMap.get(FLD_UNSPLICED_FRAGS)]))
                        .medianTpmCancer(medianCancer)
                        .percentileCancer(percentileCancer)
                        .medianTpmCohort(medianCohort)
                        .percentileCohort(percentileCohort)
                        .build());
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load gene data file({}): {}", filename, e.toString());
            return;
        }

        ISF_LOGGER.debug("sample({}) writing {} gene expression records to DB", sampleId, geneExpressions.size());
        rnaDAO.writeGeneExpressions(sampleId, geneExpressions);
    }

    private void loadNovelJunctions(final String sampleId, final IsofoxDAO rnaDAO)
    {
        if(!mConfig.loadDataType(NOVEL_JUNCTION))
            return;

        final List<NovelSpliceJunction> novelJunctions = Lists.newArrayList();

        final String sampleDataDir = mConfig.AltSjDataDir.contains("*") ?
                mConfig.AltSjDataDir.replaceAll("\\*", sampleId) : mConfig.AltSjDataDir;

        final String filename = sampleDataDir + sampleId + ISF_FILE_ID + ALT_SJ_FILE_ID;

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();
            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIMITER);

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                int geneIdIndex = fieldsIndexMap.get(FLD_GENE_ID);
                int geneName = fieldsIndexMap.get(FLD_GENE_NAME);
                int chr = fieldsIndexMap.get(FLD_CHROMOSOME);
                int posStart = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
                int posEnd = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
                int type = fieldsIndexMap.get(FLD_ALT_SJ_TYPE);
                int fragCount = fieldsIndexMap.get(FLD_ALT_SJ_FRAG_COUNT);
                int depthStart = fieldsIndexMap.get("DepthStart");
                int depthEnd = fieldsIndexMap.get("DepthEnd");
                int regionStart = fieldsIndexMap.containsKey("RegionStart") ? fieldsIndexMap.get("RegionStart") : fieldsIndexMap.get("ContextStart");
                int regionEnd = fieldsIndexMap.containsKey("RegionEnd") ? fieldsIndexMap.get("RegionEnd") : fieldsIndexMap.get("ContextEnd");
                int basesStart = fieldsIndexMap.get("BasesStart");
                int basesEnd = fieldsIndexMap.get("BasesEnd");
                int transStart = fieldsIndexMap.get("TransStart");
                int transEnd = fieldsIndexMap.get("TransEnd");

                final String geneId = items[geneIdIndex];

                if(!mConfig.processGeneId(geneId))
                    continue;

                final AltSpliceJunctionFile altSJ = AltSpliceJunctionFile.fromCsv(
                        items, geneIdIndex, geneName, chr, posStart, posEnd, type,
                        fragCount, depthStart, depthEnd, regionStart, regionEnd, basesStart, basesEnd, transStart, transEnd);

                Integer cohortFrequency = mAltSjCohortFrequency.get(altSJ.key());

                novelJunctions.add(ImmutableNovelSpliceJunction.builder()
                        .geneName(items[fieldsIndexMap.get(FLD_GENE_NAME)])
                        .chromosome(altSJ.Chromosome)
                        .junctionStart(altSJ.SpliceJunction[SE_START])
                        .junctionEnd(altSJ.SpliceJunction[SE_END])
                        .type(altSJ.Type.toString())
                        .fragmentCount(altSJ.FragmentCount)
                        .depthStart(altSJ.DepthCounts[SE_START])
                        .depthEnd(altSJ.DepthCounts[SE_END])
                        .regionStart(altSJ.RegionContexts[SE_START].toString())
                        .regionEnd(altSJ.RegionContexts[SE_END].toString())
                        .basesStart(altSJ.BaseContexts[SE_START])
                        .basesEnd(altSJ.BaseContexts[SE_END])
                        .cohortFrequency(cohortFrequency != null ? cohortFrequency : 0)
                        .build());
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load alt-SJ file({}): {}", filename, e.toString());
            return;
        }

        ISF_LOGGER.debug("sample({}) writing {} novel SJs records to DB", sampleId, novelJunctions.size());
        rnaDAO.writeNovelSpliceJunctions(sampleId, novelJunctions);
    }

    private void loadFusions(final String sampleId, final IsofoxDAO rnaDAO)
    {
        if(!mConfig.loadDataType(DataLoadType.FUSION))
            return;

        final List<RnaFusion> fusions = Lists.newArrayList();

        final String sampleDataDir = mConfig.FusionDataDir.contains("*") ?
                mConfig.FusionDataDir.replaceAll("\\*", sampleId) : mConfig.FusionDataDir;

        final String filename = sampleDataDir + sampleId + ISF_FILE_ID + PASS_FUSION_FILE_ID;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));
            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), DELIMITER);
            lines.remove(0);

            for(String line : lines)
            {
                final String[] items = line.split(DELIMITER, -1);

                String fusionName = String.format("%s_%s",
                        items[fieldsIndexMap.get(formStreamField(FLD_GENE_NAME, FS_UP))],
                        items[fieldsIndexMap.get(formStreamField(FLD_GENE_NAME, FS_DOWN))]);

                fusions.add(ImmutableRnaFusion.builder()
                        .name(fusionName)
                        .chromosomeUp(items[fieldsIndexMap.get(formStreamField(FLD_CHR, FS_UP))])
                        .chromosomeDown(items[fieldsIndexMap.get(formStreamField(FLD_CHR, FS_DOWN))])
                        .positionUp(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_POS, FS_UP))]))
                        .positionDown(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_POS, FS_DOWN))]))
                        .orientationUp(Byte.parseByte(items[fieldsIndexMap.get(formStreamField(FLD_ORIENT, FS_UP))]))
                        .orientationDown(Byte.parseByte(items[fieldsIndexMap.get(formStreamField(FLD_ORIENT, FS_DOWN))]))
                        .junctionTypeUp(items[fieldsIndexMap.get(formStreamField(FLD_JUNC_TYPE, FS_UP))])
                        .junctionTypeDown(items[fieldsIndexMap.get(formStreamField(FLD_JUNC_TYPE, FS_DOWN))])
                        .svType(items[fieldsIndexMap.get(FLD_SV_TYPE)])
                        .splitFragments(Integer.parseInt(items[fieldsIndexMap.get(FLD_SPLIT_FRAGS)]))
                        .realignedFrags(Integer.parseInt(items[fieldsIndexMap.get(FLD_REALIGN_FLAGS )]))
                        .discordantFrags(Integer.parseInt(items[fieldsIndexMap.get(FLD_DISCORD_FRAGS)]))
                        .depthUp(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_COVERAGE, FS_UP))]))
                        .depthDown(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_COVERAGE, FS_DOWN))]))
                        .maxAnchorLengthUp(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_COVERAGE, FS_UP))]))
                        .maxAnchorLengthDown(Integer.parseInt(items[fieldsIndexMap.get(formStreamField(FLD_COVERAGE, FS_DOWN))]))
                        .cohortFrequency(Integer.parseInt(items[fieldsIndexMap.get(FLD_COHORT_COUNT)]))
                        .build());
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load fusions file({}): {}", filename, e.toString());
            return;
        }

        ISF_LOGGER.debug("sample({}) writing {} fusion records to DB", sampleId, fusions.size());
        rnaDAO.writeRnaFusions(sampleId, fusions);
    }

    private void loadAltSjCohortFile()
    {
        if(!mConfig.loadDataType(NOVEL_JUNCTION))
            return;

        // isofox.driver_fusion_2212.alt_sj_cohort.csv
        // GeneId,SampleCount,Chromosome,Type,SjStart,SjEnd,

        if(mConfig.AltSjCohortFile == null || !Files.exists(Paths.get(mConfig.AltSjCohortFile)))
        {
            ISF_LOGGER.error("missing alt-SJ cohort file");
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(mConfig.AltSjCohortFile));

            String line = fileReader.readLine();
            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIMITER);

            int geneIndex = fieldsIndexMap.get(FLD_GENE_ID);
            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int sampleCountIndex = fieldsIndexMap.get("SampleCount");

            while((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER, -1);

                final String geneId = items[geneIndex];

                if(!mConfig.processGeneId(geneId))
                    continue;

                int sampleCount = Integer.parseInt(items[sampleCountIndex]);
                final String asjKey = formKey(items[chrIndex], Integer.parseInt(items[posStartIndex]), Integer.parseInt(items[posEndIndex]));

                mAltSjCohortFrequency.put(asjKey, sampleCount);
            }

            ISF_LOGGER.info("loaded alt-SJ cohort file({}) with {} sites", mConfig.AltSjCohortFile, mAltSjCohortFrequency.size());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load alt-SJ cohort file({}): {}", mConfig.AltSjCohortFile, e.toString());
            return;
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = DataLoaderConfig.createCmdLineOptions();
        final CommandLineParser parser = new DefaultParser();
        final CommandLine cmd = parser.parse(options, args);

        setLogLevel(cmd);

        IsofoxDataLoader dataLoader = new IsofoxDataLoader(cmd);

        if(!dataLoader.load())
        {
            ISF_LOGGER.info("Isofox data loading failed");
            return;
        }

        ISF_LOGGER.info("Isofox data loading complete");
    }

}
