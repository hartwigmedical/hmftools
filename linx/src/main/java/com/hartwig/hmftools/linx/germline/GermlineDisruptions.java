package com.hartwig.hmftools.linx.germline;

import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UNKNOWN;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.DOWNSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusterMetrics.findEndIndex;
import static com.hartwig.hmftools.linx.analysis.ClusterMetrics.findStartIndex;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.fusion.SvDisruptionData;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class GermlineDisruptions
{
    private final EnsemblDataCache mGeneTransCache;
    private final List<GeneData> mDriverGeneDataList;
    private final List<String> mReportableGeneIds;
    private final GermlinePonCache mGermlinePonCache;
    private final List<SvDisruptionData> mDisruptions;

    private static final int MAX_DELETE_LENGTH = 5000000;

    private static final List<ResolvedType> REPORTED_RESOLVED_TYPES = Lists.newArrayList(
            ResolvedType.DEL, ResolvedType.DUP, RECIP_INV, ResolvedType.SGL);

    public GermlineDisruptions(final LinxConfig config, final EnsemblDataCache geneTransCache, final GermlinePonCache germlinePonCache)
    {
        mGeneTransCache = geneTransCache;

        mGermlinePonCache = germlinePonCache;

        mDriverGeneDataList = config.DriverGenes.stream()
                .map(x -> geneTransCache.getGeneDataByName(x.gene()))
                .filter(x -> x != null)
                .collect(Collectors.toList());

        mReportableGeneIds = config.DriverGenes.stream()
                .filter(x -> x.reportGermlineDisruption())
                .map(x -> geneTransCache.getGeneDataByName(x.gene()))
                .filter(x -> x != null)
                .map(x -> x.GeneId)
                .collect(Collectors.toList());

        mDisruptions = Lists.newArrayList();
    }

    public final List<SvDisruptionData> getDisruptions() { return mDisruptions; }

    private boolean isReportable(final SvVarData var, final String geneId)
    {
        if(!var.getSvData().filter().equals(PASS))
            return false;

        if(!mReportableGeneIds.contains(geneId))
            return false;

        final SvCluster cluster = var.getCluster();

        if(!REPORTED_RESOLVED_TYPES.contains(cluster.getResolvedType()))
            return false;

        if(cluster.getSvCount() == 1)
        {
            if(var.type() == DEL && var.length() > MAX_DELETE_LENGTH)
                return false;
        }
        else
        {
            if(cluster.getMetrics().TotalDeleted > MAX_DELETE_LENGTH)
                return false;
        }

        return true;
    }

    public int getPonCount(final SvVarData var)
    {
        return var.getSvData().filter().equals(PON_FILTER_PON) ? mGermlinePonCache.getPonCount(var) : 0;
    }

    public void findGeneDeletions(final List<SvCluster> clusters)
    {
        mDisruptions.clear();

        for(SvCluster cluster : clusters)
        {
            if(cluster.getSvCount() == 1 && cluster.getSV(0).type() != DEL)
                continue;

            if(cluster.getSvCount() == 1)
            {
                if(cluster.getSV(0).type() != DEL)
                    continue;
            }
            else
            {
                // must be fully chained and not LINE
                if(cluster.getResolvedType() == LINE)
                    continue;

                if(!cluster.isFullyChained(true))
                    continue;
            }

            for(final Map.Entry<String, List<SvBreakend>> entry : cluster.getChrBreakendMap().entrySet())
            {
                final String chromosome = entry.getKey();
                final List<SvBreakend> breakendList = entry.getValue();

                int startIndex = findStartIndex(breakendList);
                int endIndex = findEndIndex(breakendList);

                // find stand-alone DELs and clustered deletion bridges, then look within them for driver genes which have been deleted
                for(int i = startIndex; i <= endIndex - 1; ++i)
                {
                    final SvBreakend breakend = breakendList.get(i);
                    final SvVarData var = breakend.getSV();
                    final SvBreakend nextBreakend = breakendList.get(i + 1);

                    boolean isDB = breakend.getDBLink() != null && breakend.getDBLink() == nextBreakend.getDBLink();

                    boolean isSimpleDel = !isDB && var.type() == DEL
                            && breakend.orientation() == 1 && nextBreakend.getSV() == var;

                    if(isDB || isSimpleDel)
                    {
                        int delStart = breakend.position();
                        int delEnd = nextBreakend.position();

                        if(delEnd - delStart > MAX_DELETE_LENGTH) // require a plausible deletion length
                            continue;

                        for(GeneData geneData : mDriverGeneDataList)
                        {
                            if(!geneData.Chromosome.equals(chromosome))
                                continue;

                            if(positionsWithin(geneData.GeneStart, geneData.GeneEnd, delStart, delEnd))
                            {
                                TranscriptData canonicalTrans = mGeneTransCache.getTranscriptData(geneData.GeneId, "");

                                if(canonicalTrans == null)
                                {
                                    LNX_LOGGER.error("gene({}:{}) missing canonical transcript", geneData.GeneId, geneData.GeneName);
                                    continue;
                                }

                                SvDisruptionData upDisruptionData = new SvDisruptionData(
                                        breakend.getSV(), breakend.usesStart(), geneData, canonicalTrans,
                                        new int[] { 1, canonicalTrans.exons().size() + 1 }, UNKNOWN, UPSTREAM, 1.0);

                                mDisruptions.add(upDisruptionData);

                                SvDisruptionData downDisruptionData = new SvDisruptionData(
                                        nextBreakend.getSV(), nextBreakend.usesStart(), geneData, canonicalTrans,
                                        new int[] { 1, canonicalTrans.exons().size() + 1 }, UNKNOWN, DOWNSTREAM, 1.0);

                                mDisruptions.add(downDisruptionData);
                            }
                        }
                    }
                }
            }
        }
    }

    public static String csvHeader()
    {
        return (",Filter,QualScore,PonCount,NormRefFragsStart,NormRefFragsEnd,NormVarFrags");
    }

    public List<String> formCohortData(final String sampleId, final List<SvDisruptionData> standardDisruptions)
    {
        List<String> outputLines = Lists.newArrayList();

        final List<SvDisruptionData> allDisruptions = Lists.newArrayList(mDisruptions);
        allDisruptions.addAll(standardDisruptions);

        for(final SvDisruptionData disruptionData : allDisruptions)
        {
            final SvVarData var = disruptionData.Var;
            final GeneData gene = disruptionData.Gene;

            // reassessed with specific germline rules
            disruptionData.setReportable(isReportable(var, gene.GeneId));

            StringBuilder sb = new StringBuilder();

            sb.append(String.format("%s,%s,%s,%.2f,%d",
                    sampleId, disruptionData.asCsv(), var.getSvData().filter(), var.getSvData().qualityScore(), getPonCount(var)));

            sb.append(String.format(",%d,%d,%d",
                    var.getSvData().startTumorReferenceFragmentCount(), var.getSvData().endTumorReferenceFragmentCount(),
                    var.getSvData().startTumorVariantFragmentCount()));

            outputLines.add(sb.toString());
        }

        return outputLines;
    }

}
