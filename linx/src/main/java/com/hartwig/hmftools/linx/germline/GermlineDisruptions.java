package com.hartwig.hmftools.linx.germline;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_LINX_GERMLINE;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UNKNOWN;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.DOWNSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PASS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusterMetrics.findEndIndex;
import static com.hartwig.hmftools.linx.analysis.ClusterMetrics.findStartIndex;
import static com.hartwig.hmftools.linx.annotators.PseudoGeneFinder.isPseudogeneDeletion;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_TRANS;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxGermlineSv;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.fusion.SvDisruptionData;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class GermlineDisruptions
{
    private final EnsemblDataCache mGeneTransCache;
    private final List<GeneData> mDriverGeneDataList;
    private final List<String> mReportableGeneIds;
    private final List<SvDisruptionData> mDisruptions;

    private final Set<SvVarData> mReportableSgls;

    private static final int MAX_DELETE_LENGTH = 3000000;
    private static final int MAX_SGL_MAPPED_LENGTH = 500000;
    private static final String FILTER_PSEUDOGENE = "PSEUDOGENE";

    private static final List<ResolvedType> REPORTED_RESOLVED_TYPES = Lists.newArrayList(
            ResolvedType.DEL, ResolvedType.DUP, RECIP_INV, RECIP_TRANS);

    public GermlineDisruptions(final LinxConfig config, final EnsemblDataCache geneTransCache)
    {
        mGeneTransCache = geneTransCache;

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

        mReportableSgls = Sets.newHashSet();
    }

    private boolean isReportable(final SvDisruptionData disruptionData)
    {
        final SvVarData var = disruptionData.Var;

        if(!var.getSvData().filter().equals(PASS))
            return false;

        if(!mReportableGeneIds.contains(disruptionData.Gene.GeneId))
            return false;

        final SvCluster cluster = var.getCluster();

        if(!mReportableSgls.contains(var) && !REPORTED_RESOLVED_TYPES.contains(cluster.getResolvedType()))
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

        if(disruptionData.isPseudogeneDeletion())
            return false;

        return true;
    }

    private int getPonCount(final SvVarData var)
    {
        return var.getSvData().filter().equals(PON_FILTER_PON) ? var.getSvData().ponCount() : 0;
    }

    public void findGeneDeletions(final List<SvCluster> clusters)
    {
        mDisruptions.clear();
        mReportableSgls.clear();

        for(SvCluster cluster : clusters)
        {
            if(cluster.getTypeCount(SGL) > 0)
            {
                cluster.getSVs().stream().filter(x -> x.type() == SGL).forEach(x -> checkSglMappings(x));
            }

            if(cluster.getSvCount() == 1 && cluster.getSV(0).type() != DEL)
                continue;

            checkClusterGeneDeletions(cluster);
        }
    }

    private void checkClusterGeneDeletions(final SvCluster cluster)
    {
        // must be fully chained (if not a single DEL) and not LINE
        if(cluster.getResolvedType() == LINE)
            return;

        if(cluster.getSvCount() > 1 && !cluster.isFullyChained(true))
            return;

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
                        return;

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

    private void checkSglMappings(final SvVarData var)
    {
        // look for SGLs which have mappings so as to make them a DEL or DUP candidate
        SvBreakend breakendStart = var.getBreakend(true);

        for(SglMapping mapping : var.getSglMappings())
        {
            if(!mapping.Chromosome.equals(breakendStart.chromosome()))
                continue;

            if(mapping.Orientation == breakendStart.orientation())
                continue;

            int posStart = 0;
            int posEnd = 0;
            StructuralVariantType impliedType = null;

            if(mapping.Position < breakendStart.position())
            {
                posStart = mapping.Position;
                posEnd = breakendStart.position();
                impliedType = mapping.Orientation == POS_ORIENT ? DEL : DUP;
            }
            else
            {
                posStart = breakendStart.position();
                posEnd = mapping.Position;
                impliedType = breakendStart.orientation() == POS_ORIENT ? DEL : DUP;
            }

            if(abs(mapping.Position - breakendStart.position()) > MAX_SGL_MAPPED_LENGTH)
                continue;

            for(GeneData geneData : mDriverGeneDataList)
            {
                if(!geneData.Chromosome.equals(breakendStart.chromosome()))
                    continue;

                TranscriptData canonicalTrans = mGeneTransCache.getTranscriptData(geneData.GeneId, "");

                if(canonicalTrans == null)
                {
                    LNX_LOGGER.error("gene({}:{}) missing canonical transcript", geneData.GeneId, geneData.GeneName);
                    continue;
                }

                // first check whole gene deletion
                boolean isDisruptive = false;

                if(impliedType == DEL && positionsWithin(geneData.GeneStart, geneData.GeneEnd, posStart, posEnd))
                {
                    isDisruptive = true;
                }
                else
                {
                    if(!positionWithin(posStart, canonicalTrans.TransStart, canonicalTrans.TransEnd)
                    && !positionWithin(posEnd, canonicalTrans.TransStart, canonicalTrans.TransEnd))
                    {
                        continue;
                    }

                    int exonRankStart= -1;
                    int exonRankEnd = -1;

                    for(int i = 0; i < canonicalTrans.exons().size(); ++i)
                    {
                        ExonData exon = canonicalTrans.exons().get(i);

                        if(positionWithin(posStart, exon.Start, exon.End) || positionWithin(posEnd, exon.Start, exon.End))
                        {
                            isDisruptive = true;
                            break;
                        }

                        ExonData nextExon = i < canonicalTrans.exons().size() - 1 ? canonicalTrans.exons().get(i + 1) : null;

                        if(posStart > exon.End && nextExon != null && posStart < nextExon.Start)
                        {
                            exonRankStart = min(exon.Rank, nextExon.Rank);
                        }

                        if(posEnd > exon.End && nextExon != null && posEnd < nextExon.Start)
                        {
                            exonRankEnd = min(exon.Rank, nextExon.Rank);
                        }
                    }

                    if(exonRankStart != exonRankEnd)
                    {
                        if(impliedType == DUP)
                        {
                            // cannot just duplicate the start of a gene
                            if(canonicalTrans.posStrand() && exonRankStart >= 1)
                                isDisruptive = true;
                            else if(!canonicalTrans.posStrand() && exonRankEnd >= 1)
                                isDisruptive = true;
                        }
                        else
                        {
                            isDisruptive = true;
                        }
                    }
                }

                if(isDisruptive)
                {
                    SvDisruptionData disruptionData = new SvDisruptionData(
                            breakendStart.getSV(), breakendStart.usesStart(), geneData, canonicalTrans,
                            new int[] { 1, canonicalTrans.exons().size() + 1 }, UNKNOWN, UPSTREAM, 1.0);

                    if(impliedType == DEL && isPseudogeneDeletion(var, posStart, posEnd, canonicalTrans))
                    {
                        disruptionData.markPseudogeneDeletion();
                    }

                    mDisruptions.add(disruptionData);
                    mReportableSgls.add(var);
                }
            }
        }
    }

    public void writeGermlineSVs(
            final List<SvDisruptionData> standardDisruptions, final String sampleId, final String outputDir, final DatabaseAccess dbAccess)
    {
        List<LinxGermlineSv> germlineSVs = Lists.newArrayList();
        List<DriverCatalog> drivers = Lists.newArrayList();

        populateGermlineSVs(standardDisruptions, germlineSVs, drivers);

        if(outputDir != null)
        {
            try
            {
                // write flat files for database loading
                LinxGermlineSv.write(LinxGermlineSv.generateFilename(outputDir, sampleId), germlineSVs);

                DriverCatalogFile.write(LinxDriver.generateCatalogFilename(outputDir, sampleId, false), drivers);
            }
            catch(IOException e)
            {
                LNX_LOGGER.error("failed to write germline SV file: {}", e.toString());
            }
        }

        if(dbAccess != null)
        {
            LNX_LOGGER.info("uploading {} germline SVs to database", germlineSVs.size());
            dbAccess.writeGermlineSVs(sampleId, germlineSVs);

            dbAccess.writeLinxDriverCatalog(sampleId, drivers, DRIVERS_LINX_GERMLINE);
        }
    }

    public void populateGermlineSVs(
            final List<SvDisruptionData> standardDisruptions, final List<LinxGermlineSv> germlineSVs, final List<DriverCatalog> drivers)
    {
        final List<SvDisruptionData> allDisruptions = Lists.newArrayList(mDisruptions);
        allDisruptions.addAll(standardDisruptions);

        Map<SvVarData,Set<String>> processedSvGenes = Maps.newHashMap();

        for(final SvDisruptionData disruptionData : allDisruptions)
        {
            final SvVarData var = disruptionData.Var;

            final GeneData gene = disruptionData.Gene;

            Set<String> processedGenes = processedSvGenes.get(var);

            if(processedGenes == null)
            {
                processedGenes = Sets.newHashSet();
                processedSvGenes.put(var, processedGenes);
            }
            else if(processedGenes.contains(gene.GeneName))
            {
                continue;
            }

            processedGenes.add(gene.GeneName);

            StructuralVariantData svData = var.getSvData();

            SvCluster cluster = var.getCluster();

            int ponCount = getPonCount(var);
            boolean reportable = isReportable(disruptionData);

            // TODO: switch tumor for normal fragment counts until GRIPSS writes them both

            String filters = svData.filter();

            if(disruptionData.isPseudogeneDeletion())
            {
                if(filters.equals(PASS))
                    filters = FILTER_PSEUDOGENE;
                else
                    filters += ";" + FILTER_PSEUDOGENE;
            }

            germlineSVs.add(new LinxGermlineSv(
                    var.chromosome(true), var.chromosome(false),
                    var.position(true), var.position(false),
                    var.orientation(true), var.orientation(false),
                    gene.GeneName, var.type(), filters, svData.event(), svData.qualityScore(),
                    svData.startNormalVariantFragmentCount(), svData.startNormalReferenceFragmentCount(), svData.endNormalReferenceFragmentCount(),
                    svData.startTumorVariantFragmentCount(), svData.startTumorReferenceFragmentCount(), svData.endTumorReferenceFragmentCount(),
                    svData.insertSequence(), svData.insertSequenceAlignments(), svData.insertSequenceRepeatClass(), svData.insertSequenceRepeatType(),
                    cluster.id(), cluster.getSvCount(), cluster.getResolvedType().toString(),
                    svData.startLinkedBy(), svData.endLinkedBy(), ponCount, reportable));

            if(reportable)
            {
                drivers.add(ImmutableDriverCatalog.builder()
                        .driver(DriverType.GERMLINE_DISRUPTION)
                        .category(TSG)
                        .gene(gene.GeneName)
                        .transcript(disruptionData.Transcript.TransName)
                        .isCanonical(disruptionData.Transcript.IsCanonical)
                        .chromosome(gene.Chromosome)
                        .chromosomeBand(gene.KaryotypeBand)
                        .likelihoodMethod(LikelihoodMethod.GERMLINE)
                        .driverLikelihood(1.0)
                        .missense(0)
                        .nonsense(0)
                        .splice(0)
                        .inframe(0)
                        .frameshift(0)
                        .biallelic(false)
                        .minCopyNumber(0)
                        .maxCopyNumber(0)
                        .build());
            }
        }
    }

    public static String csvHeader()
    {
        return (",Filter,QualScore,PonCount,IsPseudogene,NormRefFragsStart,NormRefFragsEnd,NormVarFrags");
    }

    public List<String> formCohortData(final String sampleId, final List<SvDisruptionData> standardDisruptions)
    {
        List<String> outputLines = Lists.newArrayList();

        final List<SvDisruptionData> allDisruptions = Lists.newArrayList(mDisruptions);
        allDisruptions.addAll(standardDisruptions);

        Set<SvVarData> processedSgls = Sets.newHashSet();

        for(final SvDisruptionData disruptionData : allDisruptions)
        {
            final SvVarData var = disruptionData.Var;

            if(var.isSglBreakend())
            {
                if(processedSgls.contains(var))
                    continue;

                processedSgls.add(var);
            }

            // reassessed with specific germline rules
            disruptionData.setReportable(isReportable(disruptionData));

            StringBuilder sb = new StringBuilder();

            sb.append(String.format("%s,%s,%s,%.2f,%d,%s",
                    sampleId, disruptionData.asCsv(), var.getSvData().filter(), var.getSvData().qualityScore(), getPonCount(var),
                    disruptionData.isPseudogeneDeletion()));

            sb.append(String.format(",%d,%d,%d",
                    var.getSvData().startTumorReferenceFragmentCount(), var.getSvData().endTumorReferenceFragmentCount(),
                    var.getSvData().startTumorVariantFragmentCount()));

            outputLines.add(sb.toString());
        }

        return outputLines;
    }

}
