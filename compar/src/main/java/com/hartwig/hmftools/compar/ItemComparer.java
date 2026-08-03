package com.hartwig.hmftools.compar;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.common.TruthsetValue;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public abstract class ItemComparer
{
    protected final ComparConfig mConfig;
    protected final List<FieldInfo> mFields;

    public ItemComparer(final ComparConfig config)
    {
        mConfig = config;
        mFields = Lists.newArrayList();
    }

    public abstract CategoryType category();

    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        return loadFromFile(sampleId, germlineSampleId, fileSources);
    }

    public List<ComparableItem> loadFromTruthset(final List<TruthsetValue> truthsetValues)
    {
        return Collections.emptyList();
    }

    public List<FieldInfo> fieldsList() { return mFields; }

    public List<String> displayFieldNames()
    {
        return mFields.stream().map(x -> x.name()).collect(Collectors.toList());
    }

    public boolean hasReportable() { return true; }

}
