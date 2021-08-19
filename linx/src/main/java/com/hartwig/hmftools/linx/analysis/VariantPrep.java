package com.hartwig.hmftools.linx.analysis;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.linx.fusion.FusionConstants.PRE_GENE_PROMOTOR_DISTANCE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.linx.cn.JcnCalcData;
import com.hartwig.hmftools.linx.cn.SvCNData;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvVarData;

public class VariantPrep
{
    public static void setSvCopyNumberData(List<SvVarData> svList, final Map<Integer, JcnCalcData> svJcnCalcDataMap,
            final Map<Integer, SvCNData[]> svIdCnDataMap, final Map<String,List<SvCNData>> chrCnDataMap)
    {
        if((svJcnCalcDataMap == null || svJcnCalcDataMap.isEmpty()) && svIdCnDataMap.isEmpty())
            return;

        List<SvCNData> cnDataList = null;
        String currentChromosome = "";
        for(final SvVarData var : svList)
        {
            if(svJcnCalcDataMap != null)
            {
                final JcnCalcData jcnData = svJcnCalcDataMap.get(var.id());
                if (jcnData != null)
                {
                    double estJcn = jcnData.JcnEstimate;
                    double estUncertainty = jcnData.JcnUncertainty;
                    var.setJcnRecalcData(estJcn - estUncertainty, estJcn + estUncertainty);
                }
            }

            final SvCNData[] cnDataPair = svIdCnDataMap.get(var.id());

            if(cnDataPair == null)
                continue;

            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(var.isSglBreakend() && be == SE_END)
                    continue;

                boolean isStart = isStart(be);

                if(!currentChromosome.equals(var.chromosome(isStart)))
                {
                    currentChromosome = var.chromosome(isStart);
                    cnDataList = chrCnDataMap.get(currentChromosome);
                }

                SvCNData cnDataPost = cnDataPair[be];

                if(cnDataList == null || cnDataPost == null)
                    continue;

                SvCNData cnDataPrev = cnDataList.get(cnDataPost.getIndex() - 1);

                var.setCopyNumberData(isStart, cnDataPrev, cnDataPost);
            }
        }
    }

    public static void setSvGeneData(
            final List<SvVarData> svList, final EnsemblDataCache ensemblDataCache, boolean applyPromotorDistance, boolean loadBreakendGenes)
    {
        int upstreamDistance = applyPromotorDistance ? PRE_GENE_PROMOTOR_DISTANCE : 0;

        if (loadBreakendGenes)
        {
            // only load transcript info for the genes covered
            final List<String> restrictedGeneIds = Lists.newArrayList();

            for (final SvVarData var : svList)
            {
                for (int be = SE_START; be <= SE_END; ++be)
                {
                    if (be == SE_END && var.isSglBreakend())
                    {
                        // special case of looking for mappings to locations containing genes so hotspot fusions can be found
                        for(final SglMapping mapping : var.getSglMappings())
                        {
                            ensemblDataCache.populateGeneIdList(restrictedGeneIds, mapping.Chromosome, mapping.Position, upstreamDistance);
                        }
                    }
                    else
                    {
                        boolean isStart = isStart(be);
                        ensemblDataCache.populateGeneIdList(restrictedGeneIds, var.chromosome(isStart), var.position(isStart), upstreamDistance);
                    }
                }
            }

            ensemblDataCache.getAlternativeGeneData().stream().filter(x -> !restrictedGeneIds.contains(x.GeneId))
                    .forEach(x -> restrictedGeneIds.add(x.GeneId));

            ensemblDataCache.loadTranscriptData(restrictedGeneIds);
        }

        // associate breakends with transcripts
        for (final SvVarData var : svList)
        {
            for (int be = SE_START; be <= SE_END; ++be)
            {
                boolean isStart = isStart(be);
                final List<BreakendGeneData> genesList = var.getGenesList(isStart);

                if (be == SE_END && var.isSglBreakend())
                {
                    // special case of looking for mappings to locations containing genes so hotspot fusions can be found
                    for(final SglMapping mapping : var.getSglMappings())
                    {
                        final List<BreakendGeneData> mappingGenes = ensemblDataCache.findGeneAnnotationsBySv(
                                var.id(), isStart, mapping.Chromosome, mapping.Position, mapping.Orientation, upstreamDistance);

                        mappingGenes.forEach(x -> x.setType(var.type()));

                        genesList.addAll(mappingGenes);
                    }
                }
                else
                {
                    genesList.addAll(ensemblDataCache.findGeneAnnotationsBySv(
                            var.id(), isStart, var.chromosome(isStart), var.position(isStart), var.orientation(isStart), upstreamDistance));

                    for (BreakendGeneData gene : genesList)
                    {
                        gene.setSvData(var.getSvData(), var.jcn());
                    }
                }
            }
        }
    }

}
