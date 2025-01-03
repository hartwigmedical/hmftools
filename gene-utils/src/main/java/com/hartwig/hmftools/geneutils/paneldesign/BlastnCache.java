package com.hartwig.hmftools.geneutils.paneldesign;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Maps;

public class BlastnCache
{
    private final String mCacheFile;
    private final Map<String,BlastnResult> mSequenceDataMap;
    private final boolean mEnabled;
    private final boolean mHasExistingData;
    private boolean mHasNewData;

    public BlastnCache(final String cacheFile)
    {
        mSequenceDataMap = Maps.newHashMap();
        mCacheFile = cacheFile;

        if(cacheFile != null)
        {
            mEnabled = true;
            loadCache(cacheFile);
            mHasExistingData = !mSequenceDataMap.isEmpty();
        }
        else
        {
            mHasExistingData = false;
            mEnabled = false;
        }

        mHasNewData = false;
    }

    public boolean enabled() { return mEnabled; }

    public void addCacheEntry(final BlastnResult result)
    {
        mSequenceDataMap.put(result.Sequence, new BlastnResult(result.Sequence, result.SumBitScore, result.ResultCount));
        mHasNewData = true;
    }

    public BlastnResult findResult(final String sequence) { return mSequenceDataMap.get(sequence); }

    private enum BlastDataColumns
    {
        Sequence,
        BitScore,
        ResultCount;
    }

    public void writeCache()
    {
        if(!mEnabled || !mHasNewData)
            return;

        try
        {
            String outputFile = mCacheFile;

            if(mHasExistingData)
            {
                outputFile = mCacheFile.replaceAll(TSV_EXTENSION, ".new" + TSV_EXTENSION);
            }

            BufferedWriter writer = createBufferedWriter(outputFile);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(BlastDataColumns.Sequence.toString());
            sj.add(BlastDataColumns.BitScore.toString());
            sj.add(BlastDataColumns.ResultCount.toString());
            writer.write(sj.toString());
            writer.newLine();

            for(BlastnResult result : mSequenceDataMap.values())
            {
                writer.write(format("%s\t%.3f\t%d", result.Sequence, result.SumBitScore, result.ResultCount));
                writer.newLine();
            }

            writer.close();
        }
        catch(IOException e)
        {
            GU_LOGGER.error("failed to write BlastN cache: {}", e.toString());
        }
    }

    public void loadCache(final String cacheFile)
    {
        if(!Files.exists(Paths.get(cacheFile)))
            return;

        try
        {
            BufferedReader reader = createBufferedReader(cacheFile);

            String header = reader.readLine();

            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, TSV_DELIM);

            int seqIndex = BlastDataColumns.Sequence.ordinal();
            int scoreIndex = BlastDataColumns.BitScore.ordinal();
            int resultCountIndex = BlastDataColumns.ResultCount.ordinal();

            String line = null;

            while((line = reader.readLine()) != null)
            {
                String[] values = line.split(TSV_DELIM, -1);

                String sequence = values[seqIndex];
                double bitScore = Double.parseDouble(values[scoreIndex]);
                int resultCount = Integer.parseInt(values[resultCountIndex]);

                mSequenceDataMap.put(sequence, new BlastnResult(sequence, bitScore, resultCount));
            }

            GU_LOGGER.info("loaded {} BlastN results from cache file {}", mSequenceDataMap.size(), cacheFile);
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to load BlastN results cache file: {}", e.toString());
            System.exit(1);
        }
    }
}
