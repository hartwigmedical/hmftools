package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.redux.DuplicateFrequency;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

public class DuplicateFreqPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final String NAME = "Duplicate read group frequency";
    private static final SourceTool SOURCE_TOOL = SourceTool.REDUX;

    private static final int MAX_DUP_READS = 100;

    public DuplicateFreqPrep(CommonPrepConfig config)
    {
        mConfig = config;
    }

    public String name() { return NAME; }
    public SourceTool sourceTool() { return SOURCE_TOOL; }

    private List<DuplicateFrequency> loadDuplicateFrequencies(String sampleId) throws IOException
    {
        String baseDir = mConfig.getReduxDir(sampleId);
        String filePath = DuplicateFrequency.generateFilename(baseDir, sampleId);
        return DuplicateFrequency.read(filePath);
    }

    private static List<Feature> normaliseAndBinCounts(List<DuplicateFrequency> dupFreqs)
    {
        long totalCount = dupFreqs.stream().mapToLong(x -> x.Count).sum();
        List<Feature> features = new ArrayList<>();
        long aboveMaxDupReadsCount = 0;

        for(DuplicateFrequency dupFreq : dupFreqs)
        {
            if(dupFreq.ReadCount < MAX_DUP_READS)
            {
                features.add(new Feature(
                        FeatureType.DUPLICATE_FREQ,
                        String.valueOf(dupFreq.ReadCount),
                        (double) dupFreq.Count / totalCount,
                        SOURCE_TOOL
                ));
            }
            else
            {
                aboveMaxDupReadsCount += dupFreq.Count;
            }
        }

        features.add(new Feature(
                FeatureType.DUPLICATE_FREQ,
                String.format("≥%s", MAX_DUP_READS),
                (double) aboveMaxDupReadsCount / totalCount,
                SOURCE_TOOL
        ));

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId) throws IOException
    {
        List<DuplicateFrequency> dupFreqs = loadDuplicateFrequencies(sampleId);
        List<Feature> features = normaliseAndBinCounts(dupFreqs);
        return features;
    }
}
