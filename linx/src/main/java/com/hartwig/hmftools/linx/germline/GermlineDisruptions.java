package com.hartwig.hmftools.linx.germline;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.driver.DriverType.GERMLINE_DISRUPTION;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.NONE;
import static com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting.VARIANT_NOT_LOST;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UNKNOWN;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UTR_3P;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.DOWNSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.linx.LinxBreakend.BREAKEND_ORIENTATION_DOWNSTREAM;
import static com.hartwig.hmftools.common.linx.LinxBreakend.BREAKEND_ORIENTATION_UPSTREAM;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.ClusterMetrics.findEndIndex;
import static com.hartwig.hmftools.linx.analysis.ClusterMetrics.findStartIndex;
import static com.hartwig.hmftools.linx.annotators.PseudoGeneFinder.isPseudogeneDeletion;
import static com.hartwig.hmftools.linx.drivers.DeletionDrivers.isHomozygousDupDisruption;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.types.ResolvedType.RECIP_INV;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.fusion.SvDisruptionData;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.linx.gene.BreakendTransData;
import com.hartwig.hmftools.linx.types.ResolvedType;
import com.hartwig.hmftools.linx.types.SglMapping;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;

public class GermlineDisruptions
{
    private final EnsemblDataCache mGeneTransCache;
    private final List<GeneData> mDriverGeneDataList;
    private final List<DriverGene> mDriverGenes;
    private final List<SvDisruptionData> mDisruptions;

    private final Set<SvVarData> mReportableSgls;

    private static final int MAX_DELETE_LENGTH = 3000000;
    private static final int MAX_PON_COUNT = 10;
    private static final int MAX_SGL_MAPPED_LENGTH = 500000;
    private static final String FILTER_PSEUDOGENE = "PSEUDOGENE";

    private static final List<ResolvedType> REPORTED_RESOLVED_TYPES = Lists.newArrayList(
            ResolvedType.DEL, ResolvedType.DUP, ResolvedType.INS, RECIP_INV);

    public GermlineDisruptions(final LinxConfig config, final EnsemblDataCache geneTransCache)
    {
        mGeneTransCache = geneTransCache;

        mDriverGenes = Lists.newArrayList();
        mDriverGeneDataList = Lists.newArrayList();

        for(DriverGene driverGene : config.DriverGenes)
        {
            if(driverGene.reportGermlineDisruption() != NONE)
                mDriverGenes.add(driverGene);

            GeneData geneData = geneTransCache.getGeneDataByName(driverGene.gene());

            if(geneData != null)
                mDriverGeneDataList.add(geneData);
        }

        mDisruptions = Lists.newArrayList();

        mReportableSgls = Sets.newHashSet();
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

                if(!isDB && !isSimpleDel)
                    return;

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
                                new int[] { 1, canonicalTrans.exons().size() + 1 }, UNKNOWN, UPSTREAM,
                                breakend.copyNumberLowSide(), breakend.copyNumber());

                        mDisruptions.add(upDisruptionData);

                        SvDisruptionData downDisruptionData = new SvDisruptionData(
                                nextBreakend.getSV(), nextBreakend.usesStart(), geneData, canonicalTrans,
                                new int[] { 1, canonicalTrans.exons().size() + 1 }, UNKNOWN, DOWNSTREAM,
                                breakend.copyNumberLowSide(), breakend.copyNumber());

                        mDisruptions.add(downDisruptionData);
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
                impliedType = mapping.Orientation == ORIENT_FWD ? DEL : DUP;
            }
            else
            {
                posStart = breakendStart.position();
                posEnd = mapping.Position;
                impliedType = breakendStart.orientation() == ORIENT_FWD ? DEL : DUP;
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
                            new int[] { 1, canonicalTrans.exons().size() + 1 }, UNKNOWN, UPSTREAM,
                            breakendStart.copyNumberLowSide(),  breakendStart.copyNumber());

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

    public void writeGermlineSVs(final List<SvDisruptionData> standardDisruptions, final String sampleId, final String outputDir)
    {
        List<LinxGermlineDisruption> germlineSVs = Lists.newArrayList();
        List<DriverCatalog> drivers = Lists.newArrayList();
        List<LinxBreakend> breakends = Lists.newArrayList();

        populateGermlineSVs(standardDisruptions, germlineSVs, breakends, drivers);

        if(outputDir != null)
        {
            try
            {
                // write flat files for database loading
                LinxGermlineDisruption.write(LinxGermlineDisruption.generateFilename(outputDir, sampleId), germlineSVs);

                LinxBreakend.write(LinxBreakend.generateFilename(outputDir, sampleId, true), breakends);

                DriverCatalogFile.write(LinxDriver.generateCatalogFilename(outputDir, sampleId, false), drivers);
            }
            catch(IOException e)
            {
                LNX_LOGGER.error("failed to write germline SV file: {}", e.toString());
            }
        }
    }

    private static void mapSvDisruptions(final List<SvDisruptionData> disruptions, final Map<SvVarData,List<SvDisruptionData>> svDisruptionsMap)
    {
        for(SvDisruptionData disruptionData : disruptions)
        {
            List<SvDisruptionData> svDisruptions = svDisruptionsMap.get(disruptionData.Var);

            if(svDisruptions == null)
            {
                svDisruptions = Lists.newArrayList();
                svDisruptionsMap.put(disruptionData.Var, svDisruptions);
            }

            svDisruptions.add(disruptionData);
        }
    }

    public void populateGermlineSVs(
            final List<SvDisruptionData> standardDisruptions, final List<LinxGermlineDisruption> germlineSVs,
            final List<LinxBreakend> breakends, final List<DriverCatalog> drivers)
    {
        // each SV may have 1 or more disruptions
        Map<SvVarData,List<SvDisruptionData>> svDisruptionsMap = Maps.newHashMap();
        mapSvDisruptions(mDisruptions, svDisruptionsMap);
        mapSvDisruptions(standardDisruptions, svDisruptionsMap);

        int breakendId = 0;

        for(Map.Entry<SvVarData,List<SvDisruptionData>> entry : svDisruptionsMap.entrySet())
        {
            SvVarData var = entry.getKey();

            StructuralVariantData svData = var.getSvData();
            SvCluster cluster = var.getCluster();

            int ponCount = var.getSvData().ponCount();

            Set<String> allFilters = Sets.newHashSet(svData.filter());
            String geneName = "";
            DriverType driverType = GERMLINE_DISRUPTION; // default if meets the driver criteria

            /*
            // check for a homozygous DUP
            if(var.type() == DUP && isHomozygousDup(var, standardDisruptions))
            {
                driverType = GERMLINE_HOM_DUP_DISRUPTION;
            }
            */

            for(SvDisruptionData disruptionData : entry.getValue())
            {
                final GeneData gene = disruptionData.Gene;
                final TranscriptData transcript = disruptionData.Transcript;

                // DELs and DUPs which straddle driver genes are not reportable
                boolean reportable = !mDisruptions.contains(disruptionData) && isReportable(disruptionData);

                if(disruptionData.isPseudogeneDeletion())
                {
                    allFilters.remove(PASS);
                    allFilters.add(FILTER_PSEUDOGENE);
                }

                boolean isUpstream = (var.orientation(disruptionData.IsStart) == ORIENT_FWD) == (gene.forwardStrand());

                ImmutableLinxBreakend.Builder builder = ImmutableLinxBreakend.builder()
                        .id(breakendId++)
                        .svId(var.id())
                        .isStart(disruptionData.IsStart)
                        .gene(gene.GeneName)
                        .geneOrientation(isUpstream ? BREAKEND_ORIENTATION_UPSTREAM : BREAKEND_ORIENTATION_DOWNSTREAM)
                        .transcriptId(transcript.TransName)
                        .canonical(transcript.IsCanonical)
                        .biotype(transcript.BioType)
                        .disruptive(false)
                        .reportedDisruption(reportable)
                        .undisruptedCopyNumber(disruptionData.UndisruptedCopyNumber)
                        .totalExonCount(transcript.exons().size());

                final BreakendGeneData breakendGene = var.getGenesList(disruptionData.IsStart).stream()
                        .filter(x -> x.geneName().equals(disruptionData.Gene.GeneName)).findFirst().orElse(null);

                if(breakendGene != null && breakendGene.canonical() != null)
                {
                    BreakendTransData breakendTransData = breakendGene.canonical();

                    builder.regionType(breakendTransData.regionType())
                        .codingType(breakendTransData.codingType())
                        .disruptive(breakendTransData.isDisruptive())
                        .exonicBasePhase(breakendTransData.Phase)
                        .nextSpliceExonRank(breakendTransData.nextSpliceExonRank())
                        .nextSpliceExonPhase(breakendTransData.Phase)
                        .nextSpliceDistance(breakendTransData.isUpstream()
                                ? breakendTransData.prevSpliceAcceptorDistance() : breakendTransData.nextSpliceAcceptorDistance())
                        .exonUp(breakendTransData.ExonUpstream)
                        .exonDown(breakendTransData.ExonDownstream);
                }
                else
                {
                    builder.regionType(disruptionData.RegionType)
                            .codingType(disruptionData.CodingType)
                            .exonicBasePhase(0)
                            .nextSpliceExonRank(disruptionData.Exons[1])
                            .nextSpliceExonPhase(0)
                            .nextSpliceDistance(0)
                            .exonUp(disruptionData.Exons[0])
                            .exonDown(disruptionData.Exons[1]);
                }

                breakends.add(builder.build());

                // add at most one driver record per gene
                if(reportable && drivers.stream().noneMatch(x -> x.gene().equals(gene.GeneName)))
                {
                    DriverGene driverGene = mDriverGenes.stream()
                            .filter(x -> x.gene().equals(disruptionData.Gene.GeneName)).findFirst().orElse(null);

                    drivers.add(ImmutableDriverCatalog.builder()
                            .driver(driverType)
                            .category(driverGene.likelihoodType())
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

                geneName = gene.GeneName;
            }

            StringJoiner filters = new StringJoiner(";");
            allFilters.forEach(x -> filters.add(x));

            germlineSVs.add(new LinxGermlineDisruption(
                    var.id(), svData.vcfIdStart(),
                    var.chromosome(true), var.chromosome(false),
                    var.position(true), var.position(false),
                    var.orientation(true), var.orientation(false),
                    var.type(), filters.toString(), svData.event(), svData.qualityScore(),
                    svData.startHomologySequence(), svData.endHomologySequence(),
                    svData.inexactHomologyOffsetStart(), svData.inexactHomologyOffsetEnd(),
                    svData.junctionCopyNumber(), svData.adjustedStartAF(), svData.adjustedEndAF(),
                    svData.adjustedStartCopyNumber(), svData.adjustedEndCopyNumber(),
                    svData.adjustedStartCopyNumberChange(), svData.adjustedEndCopyNumberChange(),
                    svData.startNormalVariantFragmentCount(), svData.startNormalReferenceFragmentCount(), svData.endNormalReferenceFragmentCount(),
                    svData.startTumorVariantFragmentCount(), svData.startTumorReferenceFragmentCount(), svData.endTumorReferenceFragmentCount(),
                    svData.insertSequence(), svData.insertSequenceAlignments(), svData.insertSequenceRepeatClass(), svData.insertSequenceRepeatType(),
                    geneName, cluster.id(), cluster.getSvCount(), cluster.getResolvedType().toString(),
                    svData.startLinkedBy(), svData.endLinkedBy(), ponCount));
        }
    }

    private boolean isHomozygousDup(final SvVarData var, final List<SvDisruptionData> standardDisruptions)
    {
        if(var.type() != DUP)
            return false;

        SvDisruptionData startDisruptionData = standardDisruptions.stream()
                .filter(x -> x.Var == var && x.IsStart).findFirst().orElse(null);

        SvDisruptionData endDisruptionData = standardDisruptions.stream()
                .filter(x -> x.Var == var && !x.IsStart).findFirst().orElse(null);

        if(startDisruptionData == null || endDisruptionData == null || startDisruptionData.Gene != endDisruptionData.Gene)
            return false;

        return isHomozygousDupDisruption(var.getBreakend(true), var.getBreakend(false), var.getSvData().junctionCopyNumber());
    }

    private boolean isReportable(final SvDisruptionData disruptionData)
    {
        SvVarData var = disruptionData.Var;

        if(!var.getSvData().filter().equals(PASS))
            return false;

        if(var.getSvData().ponCount() > MAX_PON_COUNT)
            return false;

        boolean isDelOrDup = var.type() == DEL || var.type() == DUP;

        if(disruptionData.CodingType == UTR_3P && !isDelOrDup)
            return false;

        DriverGene driverGene = mDriverGenes.stream().filter(x -> x.gene().equals(disruptionData.Gene.GeneName)).findFirst().orElse(null);

        if(driverGene == null)
            return false;

        if(driverGene.reportGermlineDisruption() == VARIANT_NOT_LOST && var.getSvData().junctionCopyNumber() < 0.1)
            return false;

        if(disruptionData.isPseudogeneDeletion())
            return false;

        final SvCluster cluster = var.getCluster();

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

        if(mReportableSgls.contains(var))
            return true;

        if(REPORTED_RESOLVED_TYPES.contains(cluster.getResolvedType()))
            return true;

        // a clustered DEL or DUP is also reportable
        return isDelOrDup;
    }
}
