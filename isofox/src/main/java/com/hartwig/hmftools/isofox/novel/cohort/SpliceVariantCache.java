package com.hartwig.hmftools.isofox.novel.cohort;

import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POSITION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_SAMPLE_ID;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatcher.COHORT_ALT_SJ_FILE;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatcher.SOMATIC_VARIANT_FILE;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatcher.SV_BREAKEND_FILE;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatcher.WRITE_VARIANT_CACHE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;

public class SpliceVariantCache
{
    private final CohortConfig mConfig;
    private final Map<String,Integer> mCohortAltSJs;
    private Map<String,List<SpliceVariant>> mSampleSpliceVariants;
    private Map<String,Map<String,List<Integer>>> mSampleSvBreakends; // sample to chromosome to list of breakend locations

    private final boolean mWriteVariantCache;

    private BufferedWriter mSomaticWriter;
    private BufferedWriter mSvBreakandWriter;

    public SpliceVariantCache(final CohortConfig config, final ConfigBuilder configBuilder)
    {
        mConfig = config;
        mCohortAltSJs = Maps.newHashMap();
        mSampleSpliceVariants = Maps.newHashMap();
        mSampleSvBreakends = Maps.newHashMap();

        mSomaticWriter = null;
        mSvBreakandWriter = null;

        boolean hasCachedFiles = false;
        if(configBuilder.hasValue(SOMATIC_VARIANT_FILE))
        {
            loadSpliceVariants(configBuilder.getValue(SOMATIC_VARIANT_FILE));
            hasCachedFiles = Files.exists(Paths.get(configBuilder.getValue(SOMATIC_VARIANT_FILE)));
        }

        if(configBuilder.hasValue(SV_BREAKEND_FILE))
        {
            loadSvBreakends(configBuilder.getValue(SV_BREAKEND_FILE));
            hasCachedFiles |= Files.exists(Paths.get(configBuilder.getValue(SV_BREAKEND_FILE)));
        }

        mWriteVariantCache = !hasCachedFiles && configBuilder.hasValue(WRITE_VARIANT_CACHE);

        if(mWriteVariantCache)
            initialiseCacheWriters();

        if(configBuilder.hasValue(COHORT_ALT_SJ_FILE))
        {
            loadCohortAltSJs(configBuilder.getValue(COHORT_ALT_SJ_FILE));
        }
    }

    public void close()
    {
        closeBufferedWriter(mSomaticWriter);
        closeBufferedWriter(mSvBreakandWriter);
    }

    public void writeVariantCache(final String sampleId, final List<SpliceVariant> spliceVariants, final Map<String,List<Integer>> svBreakends)
    {
        if(!mWriteVariantCache)
            return;

        spliceVariants.forEach(x -> writeSomaticVariant(sampleId, x));
        writeSvBreakends(sampleId, svBreakends);
    }

    public boolean hasCachedSomaticVariants() { return !mSampleSpliceVariants.isEmpty(); }
    public boolean hasCachedSvBreakends() { return !mSampleSvBreakends.isEmpty(); }

    public final List<SpliceVariant> retrieveSomaticVariants(final String sampleId)
    {
        // get and purge
        final List<SpliceVariant> spliceVariants = mSampleSpliceVariants.get(sampleId);
        mSampleSpliceVariants.remove(sampleId);
        return spliceVariants;
    }

    public final Map<String,List<Integer>> retrieveSvBreakends(final String sampleId)
    {
        final Map<String,List<Integer>> svBreakends = mSampleSvBreakends.get(sampleId);
        mSampleSvBreakends.remove(sampleId);
        return svBreakends;
    }

    private void initialiseCacheWriters()
    {
        try
        {
            final String somVarFileName = mConfig.formCohortFilename("somatic_var_cache.csv");
            mSomaticWriter = createBufferedWriter(somVarFileName, false);
            mSomaticWriter.write("SampleId,GeneName,Chromosome,Position,Type,CodingEffect,Ref,Alt,HgvsImpact,TriNucContext,LocalPhaseSet");
            mSomaticWriter.newLine();

            final String svBreakendFileName = mConfig.formCohortFilename("sv_breakend_cache.csv");
            mSvBreakandWriter = createBufferedWriter(svBreakendFileName, false);
            mSvBreakandWriter.write("SampleId,Chromosome,Position");
            mSvBreakandWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice variant file: {}", e.toString());
        }
    }

    private void writeSomaticVariant(final String sampleId, final SpliceVariant variant)
    {
        if(mSomaticWriter == null)
            return;

        try
        {
            mSomaticWriter.write(String.format("%s,%s,%s,%d,%s,%s,%s,%s,%s,%s,%d",
                    sampleId, variant.GeneName, variant.Chromosome, variant.Position, variant.Type,
                    variant.CodingEffect, variant.Ref, variant.Alt, variant.CodingImpact,
                    variant.TriNucContext, variant.LocalPhaseSet));

            mSomaticWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice variant data: {}", e.toString());
        }
    }

    private void writeSvBreakends(final String sampleId, final Map<String, List<Integer>> svBreakends)
    {
        if(mSvBreakandWriter == null)
            return;

        try
        {
            for(Map.Entry<String,List<Integer>> chrEntry : svBreakends.entrySet())
            {
                for(Integer position : chrEntry.getValue())
                {
                    mSvBreakandWriter.write(String.format("%s,%s,%d", sampleId, chrEntry.getKey(), position));
                    mSvBreakandWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write sv breakend data: {}", e.toString());
        }
    }

    public int getCohortAltSjFrequency(final AltSpliceJunctionFile altSJ)
    {
        if(mCohortAltSJs.isEmpty())
            return 0;

        Integer cohortCount = mCohortAltSJs.get(altSJ.key());
        return cohortCount != null ? cohortCount : 0;
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
            String fileDelim = inferFileDelimiter(filename);

            String line = fileReader.readLine();

            if(line == null)
                return;

            Map<String,Integer> fieldsMap = createFieldsIndexMap(line, fileDelim);

            int sampleCountIndex = fieldsMap.get("SampleCount");
            int chromosomeIndex = fieldsMap.get(FLD_CHROMOSOME);
            int sjStartPosIndex = fieldsMap.get(FLD_ALT_SJ_POS_START);
            int sjEndPosIndex = fieldsMap.get(FLD_ALT_SJ_POS_END);

            while ((line = fileReader.readLine()) != null)
            {
                String[] values = line.split(fileDelim);

                String altSjKey = NovelSpliceJunctionFile.formKey(
                        values[chromosomeIndex], Integer.parseInt(values[sjStartPosIndex]), Integer.parseInt(values[sjEndPosIndex]));

                int sampleCount = Integer.parseInt(values[sampleCountIndex]);

                mCohortAltSJs.put(altSjKey, sampleCount);
            }

            ISF_LOGGER.info("loaded {} cohort alt-SJ records", mCohortAltSJs.size());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load cohort alt-SJ data file({}): {}", filename, e.toString());
            return;
        }
    }

    private void loadSpliceVariants(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.error("invalid splice variant file({})", filename);
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if(line == null)
                return;

            String fileDelim = inferFileDelimiter(filename);
            Map<String,Integer> fieldsMap = createFieldsIndexMap(line, fileDelim);

            int sampleIdIndex = fieldsMap.get(FLD_SAMPLE_ID);

            List<SpliceVariant> spliceVariants = null;
            String currentSampleId = "";
            int variantCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                String[] values = line.split(fileDelim);
                String sampleId = values[sampleIdIndex];

                if(!sampleId.equals(currentSampleId))
                {
                    currentSampleId = sampleId;
                    spliceVariants = Lists.newArrayList();
                    mSampleSpliceVariants.put(sampleId, spliceVariants);
                }

                spliceVariants.add(SpliceVariant.fromCsv(values, fieldsMap));
                ++variantCount;
            }

            ISF_LOGGER.info("loaded {} splice variants for {} samples", variantCount, mSampleSpliceVariants.size());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load splice variant data file({}): {}", filename, e.toString());
            return;
        }
    }

    private void loadSvBreakends(final String filename)
    {
        if(!Files.exists(Paths.get(filename)))
        {
            ISF_LOGGER.error("invalid SV breakend file({})", filename);
            return;
        }

        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();

            if(line == null)
                return;

            String fileDelim = inferFileDelimiter(filename);

            Map<String,Integer> fieldsMap = createFieldsIndexMap(line, fileDelim);

            int sampleIdIndex = fieldsMap.get(FLD_SAMPLE_ID);
            int chromosomeIndex = fieldsMap.get(FLD_CHROMOSOME);
            int positionIndex = fieldsMap.get(FLD_POSITION);

            String currentSampleId = "";
            Map<String,List<Integer>> svBreakends = null;
            int breakendCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] items = line.split(fileDelim);
                String sampleId = items[sampleIdIndex];
                String chromosome = items[chromosomeIndex];
                int position = Integer.parseInt(items[positionIndex]);

                if(!sampleId.equals(currentSampleId))
                {
                    currentSampleId = sampleId;
                    svBreakends = Maps.newHashMap();
                    mSampleSvBreakends.put(sampleId, svBreakends);
                }

                List<Integer> positions = svBreakends.get(chromosome);

                if(positions == null)
                    svBreakends.put(chromosome, Lists.newArrayList(position));
                else
                    positions.add(position);

                ++breakendCount;
            }

            ISF_LOGGER.info("loaded {} SV breakeands for {} samples", breakendCount, mSampleSvBreakends.size());
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load SV breakend data file({}): {}", filename, e.toString());
            return;
        }
    }

}
