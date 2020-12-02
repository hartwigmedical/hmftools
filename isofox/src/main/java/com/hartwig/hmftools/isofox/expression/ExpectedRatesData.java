package com.hartwig.hmftools.isofox.expression;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.utils.Matrix;

public class ExpectedRatesData
{
    public final String Id;

    // equivalent of buckets - 0-N transcripts and unspliced genes
    public final List<String> Categories;

    // akin to signature names - all transcriptIds and a GeneId for each gene's unspliced region
    public final List<String> TranscriptIds;

    private Matrix mTranscriptDefinitions;

    public ExpectedRatesData(final String id)
    {
        Id = id;
        Categories = Lists.newArrayList();
        TranscriptIds = Lists.newArrayList();
        mTranscriptDefinitions = null;
    }

    public Matrix getTranscriptDefinitions() { return mTranscriptDefinitions; }

    public boolean validData()
    {
        if(Categories.isEmpty() || mTranscriptDefinitions == null)
            return false;

        if(mTranscriptDefinitions.Cols != TranscriptIds.size())
            return false;

        if(mTranscriptDefinitions.Rows != Categories.size())
            return false;

        return true;
    }

    public void initialiseTranscriptDefinitions()
    {
        if(Categories.isEmpty() || TranscriptIds.isEmpty())
            return;

        mTranscriptDefinitions = new Matrix(Categories.size(), TranscriptIds.size());
    }

    public int getTranscriptIndex(final String trans)
    {
        for(int i = 0; i < TranscriptIds.size(); ++i)
        {
            if(TranscriptIds.get(i).equals(trans))
                return i;
        }

        return -1;
    }

    public int getCategoryIndex(final String category)
    {
        for(int i = 0; i < Categories.size(); ++i)
        {
            if(Categories.get(i).equals(category))
                return i;
        }

        return -1;
    }

    public void addCategory(final String category)
    {
        if(!Categories.contains(category))
            Categories.add(category);
    }
}
