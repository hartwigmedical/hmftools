package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.appendStrList;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.EXON_RANK;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.TRANS_NAME;
import static com.hartwig.hmftools.svtools.rna_expression.RegionReadData.TRANS_REF_DELIM;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;

public class RetainedIntron
{
    public final EnsemblGeneData GeneData;

    private final List<RegionReadData> mRegions;
    private final boolean mIsStart;
    private final long mPosition;

    private int mFragmentCount;
    private int mReadDepth;

    private final String mTranscriptInfo;

    public RetainedIntron(final EnsemblGeneData geneData, final List<RegionReadData> regions, boolean isStart)
    {
        GeneData = geneData;
        mRegions = regions;
        mIsStart = isStart;

        mFragmentCount = 0;
        mReadDepth = 0;

        mPosition = isStart ? regions.get(0).start() : regions.get(0).end();

        List<String> trancriptInfo = Lists.newArrayList();

        for(RegionReadData region : mRegions)
        {
            for (final String transRef : region.getRefRegions())
            {
                String[] items = transRef.split(TRANS_REF_DELIM);
                trancriptInfo.add(String.format("%s-%s", items[TRANS_NAME], items[EXON_RANK]));
            }
        }

        mTranscriptInfo = appendStrList(trancriptInfo, ';');
    }

    public final List<RegionReadData> regions() { return mRegions; }
    public boolean isStart() { return mIsStart; }
    public long position() { return mPosition; }

    public String type() { return (GeneData.forwardStrand() == mIsStart) ? "ACCEPTOR" : "DONOR"; }

    public String transcriptInfo() { return mTranscriptInfo; }

    public void addFragmentCount() { ++mFragmentCount; }
    public int getFragmentCount() { return mFragmentCount; }

    public void addReadDepth() { ++mReadDepth; }
    public int getDepth() { return mReadDepth; }

}
