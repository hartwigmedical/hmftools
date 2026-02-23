package com.hartwig.hmftools.qsee.prep.category;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.redux.DuplicateFrequency;
import com.hartwig.hmftools.qsee.common.SampleType;
import com.hartwig.hmftools.qsee.feature.Feature;
import com.hartwig.hmftools.qsee.feature.FeatureKey;
import com.hartwig.hmftools.qsee.feature.FeatureType;
import com.hartwig.hmftools.qsee.common.MultiFieldStringBuilder;
import com.hartwig.hmftools.qsee.feature.SourceTool;
import com.hartwig.hmftools.qsee.prep.CategoryPrep;
import com.hartwig.hmftools.qsee.prep.CommonPrepConfig;

public class DuplicateFreqPrep implements CategoryPrep
{
    private final CommonPrepConfig mConfig;

    private static final SourceTool SOURCE_TOOL = SourceTool.REDUX;

    private static final String FIELD_READ_COUNT = "ReadCount";

    private static final int MAX_DUP_READS = 100;

    public DuplicateFreqPrep(CommonPrepConfig config)
    {
        mConfig = config;
    }

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
                String featureName = MultiFieldStringBuilder.formSingleField(FIELD_READ_COUNT, String.valueOf(dupFreq.ReadCount));
                FeatureKey key = new FeatureKey(featureName, FeatureType.DUPLICATE_FREQ, SOURCE_TOOL);
                features.add(new Feature(key, (double) dupFreq.Count / totalCount));
            }
            else
            {
                aboveMaxDupReadsCount += dupFreq.Count;
            }
        }

        String aboveMaxDupName = MultiFieldStringBuilder.formSingleField(FIELD_READ_COUNT, String.format("≥%s",MAX_DUP_READS));
        FeatureKey aboveMaxDupKey = new FeatureKey(aboveMaxDupName, FeatureType.DUPLICATE_FREQ, SOURCE_TOOL);
        features.add(new Feature(aboveMaxDupKey, (double) aboveMaxDupReadsCount / totalCount));

        return features;
    }

    @Override
    public List<Feature> extractSampleData(String sampleId, SampleType sampleType) throws IOException
    {
        if(sampleType != SampleType.TUMOR)
        {
            return List.of();
        }

        List<DuplicateFrequency> dupFreqs = loadDuplicateFrequencies(sampleId);
        List<Feature> features = normaliseAndBinCounts(dupFreqs);
        return features;
    }
}
