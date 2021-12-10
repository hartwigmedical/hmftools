package com.hartwig.hmftools.isofox.unmapped;

import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

import org.apache.commons.compress.utils.Lists;

public class BlatMatcher
{
    private final Map<String, List<BlatResult>> mSampleBlatResults;

    public BlatMatcher(final String filename)
    {
        mSampleBlatResults = Maps.newHashMap();

        if(filename != null)
            loadBlatResults(filename);
    }

    public List<BlatResult> getSampleBlatResults(final String sampleId) { return mSampleBlatResults.get(sampleId); }

    public BlatResult findBlatSequenceMatch(final List<BlatResult> results, final int sequenceId)
    {
        if(results == null)
            return null;

        for(BlatResult result : results)
        {
            if(result.QName.endsWith(String.format("_%d", sequenceId)))
                return result;
        }

        return null;
    }

    private void loadBlatResults(final String filename)
    {
        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));

            if(lines.size() < 6)
                return;

            // skip header
            int resultsCount = 0;
            for(int i = 5; i < lines.size(); ++i)
            {
                try
                {
                    BlatResult result = new BlatResult(lines.get(i));
                    String[] searchKey = result.QName.split("_");
                    String sampleId = searchKey[0];

                    List<BlatResult> sampleResults = mSampleBlatResults.get(sampleId);

                    if(sampleResults == null)
                    {
                        sampleResults = Lists.newArrayList();
                        mSampleBlatResults.put(sampleId, sampleResults);
                    }

                    sampleResults.add(result);
                    ++resultsCount;
                }
                catch(Exception e)
                {
                    ISF_LOGGER.error("failed to parse blat result at line({}): {}", i, e.toString());
                }
            }

            ISF_LOGGER.info("loaded {} blat results from file({})", resultsCount, filename);
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to load blat results from file({}): {}", filename, e.toString());
        }
    }



}
