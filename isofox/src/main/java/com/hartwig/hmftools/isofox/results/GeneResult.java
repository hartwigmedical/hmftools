package com.hartwig.hmftools.isofox.results;

import static com.hartwig.hmftools.isofox.common.GeneMatchType.ALT;
import static com.hartwig.hmftools.isofox.common.GeneMatchType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.GeneMatchType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.GeneMatchType.READ_THROUGH;
import static com.hartwig.hmftools.isofox.common.GeneMatchType.TOTAL;
import static com.hartwig.hmftools.isofox.common.GeneMatchType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.isofox.common.GeneMatchType.UNSPLICED;
import static com.hartwig.hmftools.isofox.common.GeneMatchType.typeAsInt;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;

import org.immutables.value.Value;

@Value.Immutable
public abstract class GeneResult
{
    public abstract EnsemblGeneData geneData();
    public abstract String collectionId();
    public abstract int intronicLength();
    public abstract int transCount();
    public abstract int totalFragments();
    public abstract int supportingTrans();
    public abstract int altFragments();
    public abstract int unsplicedFragments();
    public abstract int readThroughFragments();
    public abstract int chimericFragments();
    public abstract int duplicates();
    public abstract double unsplicedAlloc();
    public abstract double fitResiduals();

    public abstract List<TranscriptResult> transcriptResults();

    public static GeneResult createGeneResults(final GeneCollection geneCollection, final GeneReadData geneReadData, final List<TranscriptResult> transResults)
    {
        final EnsemblGeneData geneData = geneReadData.GeneData;

        long geneLength = geneData.length();
        long exonicLength = geneReadData.calcExonicRegionLength();

        final int[] fragmentCounts = geneReadData.getCounts();

        double unsplicedAlloc = geneCollection.getFitAllocation(geneData.GeneId);

        return ImmutableGeneResult.builder()
                .geneData(geneData)
                .collectionId(geneCollection.chrId())
                .intronicLength((int)(geneLength - exonicLength))
                .transCount(geneReadData.getTranscripts().size())
                .totalFragments(fragmentCounts[typeAsInt(TOTAL)])
                .supportingTrans(fragmentCounts[typeAsInt(TRANS_SUPPORTING)])
                .altFragments(fragmentCounts[typeAsInt(ALT)])
                .unsplicedFragments(fragmentCounts[typeAsInt(UNSPLICED)])
                .readThroughFragments(fragmentCounts[typeAsInt(READ_THROUGH)])
                .chimericFragments(fragmentCounts[typeAsInt(CHIMERIC)])
                .duplicates(fragmentCounts[typeAsInt(DUPLICATE)])
                .unsplicedAlloc(unsplicedAlloc)
                .fitResiduals(geneReadData.getFitResiduals())
                .transcriptResults(transResults)
                .build();
    }

}
