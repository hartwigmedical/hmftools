package com.hartwig.hmftools.purple.gene;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegionUtils;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;

public final class GeneCopyNumberFactory
{
    public static List<GeneCopyNumber> geneCopyNumbers(
            final EnsemblDataCache geneTransCache,
            final List<PurpleCopyNumber> somaticCopyNumbers, final List<PurpleCopyNumber> germlineDeletions)
    {
        final List<GeneCopyNumber> result = Lists.newArrayList();

        for(List<GeneData> geneDataList : geneTransCache.getChrGeneDataMap().values())
        {
            for(GeneData geneData : geneDataList)
            {
                List<TranscriptData> transDataList = geneTransCache.getTranscripts(geneData.GeneId);

                for(TranscriptData tranData : transDataList)
                {
                    HmfTranscriptRegion region = HmfTranscriptRegionUtils.fromTranscript(geneData, tranData);

                    final GeneCopyNumberBuilder builder = new GeneCopyNumberBuilder(region);

                    RegionZipper.zip(somaticCopyNumbers, region.exons(), builder);
                    RegionZipper.zip(germlineDeletions, region.exons(), builder);

                    GeneCopyNumber geneCopyNumber = builder.build();

                    if(geneCopyNumber.totalRegions() > 0)
                        result.add(geneCopyNumber);
                }
            }
        }

        return result;
    }
}
