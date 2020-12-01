package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.sigs.DataUtils.convertList;
import static com.hartwig.hmftools.common.sigs.Percentiles.calcPercentileValues;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.SPLICE_SITE_PERCENTILES;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.expression.cohort.TransExpressionDistribution.DISTRIBUTION_SIZE;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatcher.COHORT_ALT_SJ_FILE;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.DELIMITER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionType;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;

public class SpliceSiteCache
{
    private final CohortConfig mConfig;

    private final Map<String, Set<Integer>> mCohortAltSJPositions;

    private final Map<String,Map<Integer,List<Double>>> mSpliceSiteRates;

    public SpliceSiteCache(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mCohortAltSJPositions = Maps.newHashMap();
        mSpliceSiteRates = Maps.newHashMap();

        if(cmd.hasOption(COHORT_ALT_SJ_FILE))
            loadCohortAltSJs(cmd.getOptionValue(COHORT_ALT_SJ_FILE));
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(COHORT_ALT_SJ_FILE, true, "Cohort frequency for alt SJs");
    }

    public void createPercentiles()
    {
        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, SPLICE_SITE_PERCENTILES, filenames))
            return;

        final Map<String,Integer> fieldsIndexMap = Maps.newHashMap();

        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path sampleFile = filenames.get(i);

            loadFile(sampleId, sampleFile, fieldsIndexMap);
        }

        ISF_LOGGER.info("loaded {} sample splice-site files", mConfig.SampleData.SampleIds.size());

        // write a cohort file
        writePercentiles();
    }

    public void loadFile(final String sampleId, final Path filename, final Map<String,Integer> fieldsIndexMap)
    {
        try
        {
            final List<String> lines = Files.readAllLines(filename);

            if(fieldsIndexMap.isEmpty())
                fieldsIndexMap.putAll(createFieldsIndexMap(lines.get(0), DELIMITER));

            int chromosomeIndex = fieldsIndexMap.get("Chromosome");
            int splicePosIndex = fieldsIndexMap.get("SpliceSitePosition");
            int traverseFragsIndex = fieldsIndexMap.get("TraverseFrags");
            int supportFragsIndex = fieldsIndexMap.get("SupportFrags");

            lines.remove(0);

            int spliceSiteCount = 0;

            for(final String line : lines)
            {
                // GeneSetId,Chromosome,SpliceSitePosition,TraverseFrags,SupportFrags,SkipFrags
                final String[] items = line.split(DELIMITER);
                String chromosome = items[chromosomeIndex];
                int splicePosition = Integer.parseInt(items[splicePosIndex]);

                if(!hasAltSpliceSitePosition(chromosome, splicePosition))
                    continue;

                int traverseFrags = Integer.parseInt(items[traverseFragsIndex]);
                int supportFrags = Integer.parseInt(items[supportFragsIndex]);

                if(traverseFrags == 0 && supportFrags == 0)
                    continue;

                double supportRate = supportFrags / (double)(traverseFrags + supportFrags);

                addSpliceSiteSupport(chromosome, splicePosition, supportRate);
                ++spliceSiteCount;
            }

            ISF_LOGGER.debug("sample({}) loaded {} splice-sites", sampleId, spliceSiteCount);

        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load splice site file({}): {}", filename.toString(), e.toString());
        }
    }

    private void addSpliceSiteSupport(final String chromosome, int splicePosition, double supportRate)
    {
        Map<Integer,List<Double>> positionsMap = mSpliceSiteRates.get(chromosome);

        if(positionsMap == null)
        {
            positionsMap = Maps.newHashMap();
            mSpliceSiteRates.put(chromosome, positionsMap);
        }

        List<Double> supportRates = positionsMap.get(splicePosition);

        if(supportRates == null)
        {
            positionsMap.put(splicePosition, Lists.newArrayList(supportRate));
            return;
        }

        int index = 0;

        while(index < supportRates.size())
        {
            if(supportRate < supportRates.get(index))
                break;

            ++index;
        }

        supportRates.add(index, supportRate);
    }

    private boolean processType(AltSpliceJunctionType type)
    {
        switch(type)
        {
            case SKIPPED_EXONS:
            case NOVEL_3_PRIME:
            case NOVEL_5_PRIME:
            case NOVEL_EXON:
            case MIXED_TRANS:
            case EXON_INTRON:
                return true;
            default:
                return false;
        }
    }

    private void loadCohortAltSJs(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.error("invalid cohort alt-SJ file({})", filename);
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if (line == null)
                return;

            final Map<String,Integer> fieldsMap = createFieldsIndexMap(line, DELIMITER);

            int typeIndex = fieldsMap.get("Type");
            int chromosomeIndex = fieldsMap.get("Chromosome");
            int sjStartPosIndex = fieldsMap.get("SjStart");
            int sjEndPosIndex = fieldsMap.get("SjEnd");

            int asjCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(DELIMITER);

                AltSpliceJunctionType type = AltSpliceJunctionType.valueOf(items[typeIndex]);

                if(!processType(type))
                    continue;

                String chromosome = items[chromosomeIndex];
                int asjPosStart = Integer.parseInt(items[sjStartPosIndex]);
                int asjPosEnd = Integer.parseInt(items[sjEndPosIndex]);

                Set<Integer> positions = mCohortAltSJPositions.get(chromosome);

                if(positions == null)
                {
                    positions = Sets.newHashSet();
                    mCohortAltSJPositions.put(chromosome, positions);
                }

                positions.add(asjPosStart);
                positions.add(asjPosEnd);

                ++asjCount;
            }

            ISF_LOGGER.info("loaded {} cohort alt-SJ records", asjCount);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load cohort alt-SJ data file({}): {}", filename, e.toString());
            return;
        }
    }

    private boolean hasAltSpliceSitePosition(final String chromosome, int splicePosition)
    {
        final Set<Integer> positions = mCohortAltSJPositions.get(chromosome);

        if(positions == null)
            return false;

        return positions.contains(splicePosition);
    }

    private void writePercentiles()
    {
        final double[] percentileValues = new double[DISTRIBUTION_SIZE];

        try
        {
            final String outputFileName = mConfig.formCohortFilename("splice_site_support_perc.csv");
            BufferedWriter writer = createBufferedWriter(outputFileName, false);

            writer.write("Chromosome,SplicePosition");

            for(int i = 0; i < DISTRIBUTION_SIZE; ++i)
            {
                writer.write(String.format(",Pct_%d", i));
            }

            writer.newLine();

            for(Map.Entry<String,Map<Integer,List<Double>>> chrEntry : mSpliceSiteRates.entrySet())
            {
                final String chromosome = chrEntry.getKey();

                for(Map.Entry<Integer,List<Double>> posEntry : chrEntry.getValue().entrySet())
                {
                    final int splicePosition = posEntry.getKey();

                    final double[] supportRates = convertList(posEntry.getValue());
                    calcPercentileValues(supportRates, percentileValues);

                    writer.write(String.format("%s,%d", chromosome, splicePosition));

                    for(int i = 0; i < DISTRIBUTION_SIZE; ++i)
                    {
                        writer.write(String.format(",%.3f", percentileValues[i]));
                    }

                    writer.newLine();
                }
            }

        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice site percentiles file: {}", e.toString());
        }
    }

}
