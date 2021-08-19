package com.hartwig.hmftools.linx.visualiser.file;

import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.FUSION;

import java.util.List;

import com.google.common.collect.Lists;

public class VisSampleData
{
    private String mSampleId;
    private final List<VisGeneData> mGeneData;
    private final List<VisFusionFile> mFusions;

    public VisSampleData()
    {
        mSampleId = "";
        mGeneData = Lists.newArrayList();
        mFusions = Lists.newArrayList();
    }

    public void setSampleId(final String sampleId)
    {
        mSampleId = sampleId;
        mGeneData.clear();
        mFusions.clear();
    }

    public String sampleId() { return mSampleId; }
    public final List<VisGeneData> getGeneData() { return mGeneData; }
    public final List<VisFusionFile> getFusions() { return mFusions; }

    public void addGeneExonData(int clusterId, final String geneId, final String geneName, final String transName, int transId,
            final String chromosome, final VisGeneAnnotationType annotationType)
    {
        addGeneExonData(new VisGeneData(clusterId, geneId, geneName, transName, transId, chromosome, annotationType));
    }

    public void addGeneExonData(final VisGeneData geneData)
    {
        // no duplicates for the same clusterId - favour fusions over other types where transcripts differ
        for(VisGeneData otherGeneData : mGeneData)
        {
            if(otherGeneData.GeneId.equals(geneData.GeneId) && otherGeneData.ClusterId == geneData.ClusterId)
            {
                if(otherGeneData.AnnotationType != FUSION && geneData.AnnotationType == FUSION)
                {
                    mGeneData.remove(otherGeneData);
                    mGeneData.add(geneData);
                }

                return;
            }
        }

        mGeneData.add(geneData);
    }

    public void addFusions(final List<VisFusionFile> fusions) { mFusions.addAll(fusions); }

}
