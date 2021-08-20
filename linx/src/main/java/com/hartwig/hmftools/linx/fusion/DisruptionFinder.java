package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.abs;
import static java.lang.Math.max;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.gene.TranscriptCodingType.UNKNOWN;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.DOWNSTREAM;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.UPSTREAM;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.PON_FILTER_PON;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM;
import static com.hartwig.hmftools.linx.analysis.ClusterMetrics.findEndIndex;
import static com.hartwig.hmftools.linx.analysis.ClusterMetrics.findStartIndex;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.DISRUPTION;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.fusion.BreakendGeneData;
import com.hartwig.hmftools.common.fusion.BreakendTransData;
import com.hartwig.hmftools.linx.CohortFileInterface;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.CohortDataWriter;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.germline.GermlinePonCache;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisSampleData;

public class DisruptionFinder implements CohortFileInterface
{
    private final EnsemblDataCache mGeneTransCache;
    private final List<GeneData> mDisruptionGenes;
    private final VisSampleData mVisSampleData;

    private final List<SvDisruptionData> mDisruptions;
    private final Map<BreakendTransData,String> mRemovedDisruptions; // cached for diagnostic purposes

    private final boolean mIsGermline;
    private final GermlinePonCache mGermlinePonCache;

    private final CohortDataWriter mCohortDataWriter;

    public static final int MAX_NON_DISRUPTED_CHAIN_LENGTH = 5000;

    public DisruptionFinder(
            final LinxConfig config, final EnsemblDataCache geneTransCache,
            final GermlinePonCache germlinePonCache, final CohortDataWriter cohortDataWriter, final VisSampleData visSampleData)
    {
        mGeneTransCache = geneTransCache;

        mIsGermline = config.IsGermline;
        mGermlinePonCache = germlinePonCache;

        mDisruptionGenes = disruptionGeneIds(config.DriverGenes, !mIsGermline, geneTransCache);

        mVisSampleData = visSampleData;
        mCohortDataWriter = cohortDataWriter;

        mDisruptions = Lists.newArrayList();
        mRemovedDisruptions = Maps.newHashMap();
    }

    public static List<GeneData> disruptionGeneIds(
            final List<DriverGene> driverGenes, boolean onlyReportable, final EnsemblDataCache geneTransCache)
    {
        return driverGenes.stream()
                .filter(x -> !onlyReportable || x.reportDisruption())
                .map(x -> geneTransCache.getGeneDataByName(x.gene()))
                .filter(x -> x != null)
                .collect(Collectors.toList());
    }

    public final List<SvDisruptionData> getDisruptions() { return mDisruptions; }

    public boolean matchesDisruptionGene(final BreakendGeneData gene)
    {
        return mDisruptionGenes.stream().anyMatch(x -> gene.StableId.equals(x.GeneId));
    }

    public void addDisruptionGene(final GeneData geneData)
    {
        mDisruptionGenes.add(geneData);
    }

    public void markTranscriptsDisruptive(final List<SvVarData> svList)
    {
        mRemovedDisruptions.clear();

        for(final SvVarData var : svList)
        {
            markTranscriptsDisruptive(var);

            // inferred SGLs are always non-disruptive
            if (var.isInferredSgl())
            {
                var.getGenesList(true).stream().forEach(x -> x.transcripts().stream().forEach(y -> y.setIsDisruptive(false)));
            }
        }
    }

    private static final String NON_DISRUPT_REASON_SIMPLE_SV = "SimpleSV";
    private static final String NON_DISRUPT_REASON_LINE = "LINE";

    private void markTranscriptsDisruptive(final SvVarData var)
    {
        final List<BreakendGeneData> genesStart = var.getGenesList(true);
        final List<BreakendGeneData> genesEnd = var.getGenesList(false);

        if(genesStart.isEmpty() && genesEnd.isEmpty())
            return;

        /* test the breakend to see if:
            - it isn't chained - revert to simple disrupted definitions
            - it is chained
                - the chain returns to the same intron
                - the chain traverses a splice acceptor in any direction
        */

        final SvCluster cluster = var.getCluster();

        // set the undisrupted copy number against all canonical transcripts
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(se == SE_END && var.isSglBreakend())
                continue;

            final SvBreakend breakend = var.getBreakend(se);
            double undisruptedCopyNumber = getUndisruptedCopyNumber(breakend);

            final List<BreakendGeneData> svGenes = se == SE_START ? genesStart : genesEnd;

            for(BreakendGeneData gene : svGenes)
            {
                BreakendTransData canonicalTrans = gene.canonical();

                if(canonicalTrans != null)
                    canonicalTrans.setUndisruptedCopyNumber(undisruptedCopyNumber);

                // line clusters can insert into an intron and look disruptive if a single breakend is involved,
                // but are only inserting a (non-disruptive) shard
                if(cluster.getResolvedType() == LINE)
                {
                    gene.transcripts().stream()
                            .filter(x -> x.isIntronic())
                            .filter(x -> x.isDisruptive())
                            .forEach(x -> markNonDisruptiveTranscript(x, NON_DISRUPT_REASON_LINE));
                }
            }
        }

        if(var.isSglBreakend())
            return;

        boolean isSimpleSV = var.isSimpleType();

        if(isSimpleSV)
        {
            markNonDisruptiveGeneTranscripts(var, genesStart, genesEnd, NON_DISRUPT_REASON_SIMPLE_SV);
        }

        final List<SvChain> chains = cluster.findChains(var);

        // first test each breakend in turn to see if it forms a TI wholly within an intron and where the other breakends
        // of the TI are non-genic

        for(int se = SE_START; se <= SE_END; ++se)
        {
            final SvBreakend breakend = var.getBreakend(se);

            final List<BreakendGeneData> svGenes = se == SE_START ? genesStart : genesEnd;

            List<BreakendTransData> transList = getDisruptedTranscripts(svGenes);

            if(transList.isEmpty())
                continue;

            boolean otherBreakendNonGenic = se == SE_START ? genesEnd.isEmpty() : genesStart.isEmpty();

            final List<LinkedPair> links = var.getLinkedPairs(breakend.usesStart());

            for(final LinkedPair pair : links)
            {
                final SvBreakend otherSvBreakend = pair.getOtherBreakend(breakend);

                List<BreakendGeneData> otherGenes = otherSvBreakend.getSV().getGenesList(otherSvBreakend.usesStart());
                List<BreakendTransData> otherTransList = getDisruptedTranscripts(otherGenes);

                boolean otherSvOtherBreakendNonGenic = otherSvBreakend.getSV().getGenesList(!otherSvBreakend.usesStart()).isEmpty();

                if(otherBreakendNonGenic && otherSvOtherBreakendNonGenic)
                {
                    if(markNonDisruptiveTranscripts(transList, otherTransList, "IntronicSection"))
                    {
                        LNX_LOGGER.debug("pair({}) length({}) fully intronic)", pair, pair.length());
                        removeNonDisruptedTranscripts(transList);
                    }
                }
            }

            // next test for a breakend which returns to the same intro via a chain with the same orientation
            // without passing through any splice acceptors
            if(!transList.isEmpty())
            {
                for(SvChain chain : chains)
                {
                    checkChainedTranscripts(breakend, chain, transList);
                }
            }
        }
    }

    private static List<BreakendTransData> getDisruptedTranscripts(final List<BreakendGeneData> genes)
    {
        List<BreakendTransData> transList = Lists.newArrayList();

        for (final BreakendGeneData gene : genes)
        {
            transList.addAll(gene.transcripts().stream().filter(x -> x.isDisruptive()).collect(Collectors.toList()));
        }

        return transList;
    }

    private static void removeNonDisruptedTranscripts(final List<BreakendTransData> transList)
    {
        int index = 0;

        while(index < transList.size())
        {
            if(transList.get(index).isDisruptive())
                ++index;
            else
                transList.remove(index);
        }
    }

    public static boolean isDisruptiveInTranscript(final SvVarData var, final TranscriptData transcript)
    {
        // return false if both breakends are not disruptive in this transcript
        if(var.isSglBreakend())
            return true;

        for(int se = SE_START; se <= SE_END; ++se)
        {
            BreakendTransData matchingTrans = var.getGenesList(isStart(se)).stream()
                    .filter(x -> x.StableId.equals(transcript.GeneId))
                    .map(x -> x.transcripts().stream().filter(y -> y.transId() == transcript.TransId).findFirst().orElse(null))
                    .filter(x -> x != null).findFirst().orElse(null);

            if(matchingTrans != null && matchingTrans.isDisruptive())
                return true;
        }

        return false;
    }

    public boolean isNonDisruptiveChainedPair(final BreakendTransData upTrans, final List<BreakendGeneData> downGenesList)
    {
        if(upTrans.regionType() == EXONIC)
            return false;

        // check whether these breakends may belong to the same non-disrupted segment in any upstream transcript
        final BreakendGeneData downGene = downGenesList.stream()
                .filter(x -> x.StableId.equals(upTrans.gene().StableId)).findFirst().orElse(null);

        if(downGene == null)
            return false;

        for(BreakendTransData downTrans : downGene.transcripts())
        {
            if(downTrans.transId() != upTrans.transId())
                continue;

            if(downTrans.regionType() == upTrans.regionType() && downTrans.ExonUpstream == upTrans.ExonUpstream)
                return true;
        }

        return false;
    }

    private void checkChainedTranscripts(final SvBreakend breakend, final SvChain chain, final List<BreakendTransData> transList)
    {
        // an SV whose breakends are not both within the same intron disrupts those transcripts unless
        // a) one breakend isn't genic and the other forms a TI whole within an intron OR
        // b) both breakends are in chained sections which come back to the same intro with correct orientation and
        // without traversing a splice acceptor OR

        final List<LinkedPair> links = chain.getLinkedPairs();

        for(int i = 0; i < links.size(); ++i)
        {
            int startIndex = 0;
            boolean traverseUp = false;

            if(i == 0 && chain.getOpenBreakend(true) == breakend)
            {
                startIndex = -1;
                traverseUp = true;
            }
            else if(i == links.size() - 1 && chain.getOpenBreakend(false) == breakend)
            {
                startIndex = links.size();
                traverseUp = false;
            }
            else if(links.get(i).hasBreakend(breakend))
            {
                startIndex = i;
                traverseUp = links.get(i).secondBreakend() == breakend;
            }
            else
            {
                continue;
            }

            int index = startIndex;
            long chainLength = 0;

            while(true)
            {
                index += traverseUp ? 1 : -1;

                if(index < 0 || index >= links.size())
                    break;

                final LinkedPair nextPair = links.get(index);

                // does this next link traverse another splice acceptor?
                boolean traversesGene = pairTraversesGene(nextPair, 0, false);

                if(traversesGene)
                    return;

                chainLength += nextPair.length();

                if(chainLength > MAX_NON_DISRUPTED_CHAIN_LENGTH)
                    return;

                // no need to chain this breakend any further as soon as any of its chains have crossed another splice acceptor

                // does it return to the same intron and correct orientation for any transcripts
                final SvBreakend nextBreakend = traverseUp ?
                        nextPair.secondBreakend().getOtherBreakend() : nextPair.firstBreakend().getOtherBreakend();

                if(nextBreakend == null)
                    continue;

                if(nextBreakend.orientation() == breakend.orientation() || !nextBreakend.chromosome().equals(breakend.chromosome()))
                    continue;

                List<BreakendGeneData> otherGenes = nextBreakend.getSV().getGenesList(nextBreakend.usesStart());

                if(otherGenes.isEmpty())
                    continue;

                List<BreakendTransData> otherTransList = getDisruptedTranscripts(otherGenes);

                if(!otherTransList.isEmpty())
                {
                    String contextInfo = String.format("SameIntronNoSPA;%d-%d", abs(index - startIndex), chainLength);

                    if (markNonDisruptiveTranscripts(transList, otherTransList, contextInfo))
                    {
                        removeNonDisruptedTranscripts(transList);

                        LNX_LOGGER.debug("breakends({} & {}) return to same intron, chain({}) links({}) length({})",
                                breakend, nextBreakend, chain.id(), abs(index - startIndex), chainLength);

                        if (transList.isEmpty())
                            return;
                    }
                }
            }
        }
    }

    private void markNonDisruptiveGeneTranscripts(
            final SvVarData var, final List<BreakendGeneData> genesStart, final List<BreakendGeneData> genesEnd, final String context)
    {
        for(final BreakendGeneData geneStart : genesStart)
        {
            final BreakendGeneData geneEnd = genesEnd.stream()
                    .filter(x -> x.StableId.equals(geneStart.StableId)).findFirst().orElse(null);

            if(geneEnd == null)
                continue;

            if(var.type() == DUP)
            {
                markNonDisruptiveDups(
                        geneStart.isUpstream() ? geneStart.transcripts() : geneEnd.transcripts(),
                        !geneEnd.isUpstream() ? geneEnd.transcripts() : geneStart.transcripts());
            }

            markNonDisruptiveTranscripts(geneStart.transcripts(), geneEnd.transcripts(), context);
        }
    }

    private boolean markNonDisruptiveTranscripts(final List<BreakendTransData> transList1, final List<BreakendTransData> transList2, final String context)
    {
        boolean foundMatchingTrans = false;

        for (final BreakendTransData trans1 : transList1)
        {
            final BreakendTransData trans2 = transList2.stream()
                    .filter(x -> x.transName().equals(trans1.transName())).findFirst().orElse(null);

            if(trans2 == null)
                continue;

            if(trans1.ExonUpstream == trans2.ExonUpstream && !trans1.isExonic() && !trans2.isExonic())
            {
                foundMatchingTrans = true;

                markNonDisruptiveTranscript(trans1, context);
                markNonDisruptiveTranscript(trans2, context);
            }
        }

        return foundMatchingTrans;
    }

    private boolean markNonDisruptiveDups(final List<BreakendTransData> upstreamTransList, final List<BreakendTransData> downstreamTransList)
    {
        // special case of a DUP around the 1st exon which since it doesn't have a splice acceptor does not change the transcript
        boolean foundMatchingTrans = false;

        for (final BreakendTransData upTrans : upstreamTransList)
        {
            final BreakendTransData downTrans = downstreamTransList.stream()
                    .filter(x -> x.transName().equals(upTrans.transName())).findFirst().orElse(null);

            if(downTrans == null)
                continue;

            if(upTrans.ExonUpstream == 1 && downTrans.ExonDownstream <= 2 && !upTrans.isExonic())
            {
                foundMatchingTrans = true;

                markNonDisruptiveTranscript(upTrans, NON_DISRUPT_REASON_SIMPLE_SV);
                markNonDisruptiveTranscript(downTrans, NON_DISRUPT_REASON_SIMPLE_SV);
            }
        }

        return foundMatchingTrans;
    }

    private void markNonDisruptiveTranscript(final BreakendTransData transcript, final String context)
    {
        transcript.setIsDisruptive(false);
        registerNonDisruptedTranscript(transcript, context);
    }

    private void registerNonDisruptedTranscript(final BreakendTransData transcript, final String context)
    {
        if(mIsGermline)
            return;

        if(mRemovedDisruptions.containsKey(transcript))
            return;

        if(!transcript.isCanonical() || !matchesDisruptionGene(transcript.gene()))
            return;

        LNX_LOGGER.debug("excluding gene({}) svId({}) reason({})", transcript.geneName(), transcript.gene().id(), context);

        mRemovedDisruptions.put(transcript, context);
    }

    public boolean pairTraversesGene(final LinkedPair pair, int fusionDirection, boolean isPrecodingUpstream)
    {
        // for this pair to not affect the fusion, the section it traverses cannot cross any gene's splice acceptor
        // with the same strand direction unless that is part of a fully traversed non-coding 5' exon

        int lowerPos = pair.getBreakend(true).position();
        int upperPos = pair.getBreakend(false).position();

        List<GeneData> geneDataList = mGeneTransCache.getChrGeneDataMap().get(pair.chromosome());

        if(geneDataList == null)
            return false;

        for(GeneData geneData : geneDataList)
        {
            if(lowerPos > geneData.GeneEnd)
                continue;

            if(upperPos < geneData.GeneStart)
                break;

            if(fusionDirection == 0 || geneData.Strand == fusionDirection)
            {
                // check whether a splice acceptor is encountered within this window
                List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

                if(transDataList == null)
                    continue;

                for(final TranscriptData transData : transDataList)
                {
                    for (final ExonData exonData : transData.exons())
                    {
                        if (exonData.Rank == 1)
                            continue;

                        if ((geneData.Strand == 1 && lowerPos <= exonData.Start && upperPos >= exonData.Start)
                        || (geneData.Strand == -1 && lowerPos <= exonData.End && upperPos >= exonData.End))
                        {
                            // allow an exon to be fully traversed if the upstream transcript is pre-coding
                            if (isPrecodingUpstream && lowerPos <= exonData.Start && upperPos >= exonData.End)
                            {
                                if (geneData.Strand == 1 && (transData.CodingStart == null || upperPos < transData.CodingStart))
                                    continue;
                                else if (geneData.Strand == -1 && (transData.CodingEnd == null || lowerPos > transData.CodingEnd))
                                    continue;
                            }

                            LNX_LOGGER.trace("pair({}) direction({}) traverses splice acceptor({} {}) exon(rank{} pos={})",
                                    pair.toString(), fusionDirection, geneData.GeneName, transData.TransName,
                                    exonData.Rank, exonData.Start, exonData.End);

                            return true;
                        }
                    }
                }
            }
        }

        return false;
    }

    public void findReportableDisruptions(final List<SvVarData> svList)
    {
        mDisruptions.clear();

        for (final SvVarData var : svList)
        {
            for (int be = SE_START; be <= SE_END; ++be)
            {
                if (be == SE_END && var.isSglBreakend())
                    continue;

                final List<BreakendGeneData> tsgGenesList = var.getGenesList(isStart(be)).stream()
                        .filter(x -> matchesDisruptionGene(x)).collect(Collectors.toList());

                if(tsgGenesList.isEmpty())
                    continue;

                for(BreakendGeneData gene : tsgGenesList)
                {
                    List<BreakendTransData> reportableDisruptions = gene.transcripts().stream()
                            .filter(BreakendTransData::isCanonical)
                            .filter(BreakendTransData::isDisruptive)
                            .collect(Collectors.toList());

                    for(BreakendTransData transcript : reportableDisruptions)
                    {
                        LNX_LOGGER.debug("var({}) breakend({}) gene({}) transcript({}) is disrupted, cnLowside({})",
                                var.id(), var.getBreakend(be), gene.GeneName, transcript.transName(),
                                formatJcn(transcript.undisruptedCopyNumber()));

                        transcript.setReportableDisruption(true);

                        SvDisruptionData disruptionData = new SvDisruptionData(
                                var,  gene.isStart(), true, gene.getGeneData(), transcript.TransData,
                                new int[] { transcript.ExonUpstream, transcript.ExonDownstream },
                                transcript.codingType(), transcript.regionType(), transcript.undisruptedCopyNumber());

                        mDisruptions.add(disruptionData);

                        if(mVisSampleData != null)
                        {
                            mVisSampleData.addGeneExonData(
                                    var.getCluster().id(), gene.StableId, gene.GeneName,
                                    "", 0, gene.chromosome(), DISRUPTION);
                        }
                    }
                }
            }
        }
    }

    private static double getUndisruptedCopyNumber(final SvBreakend breakend)
    {
        double cnLowSide = breakend.copyNumberLowSide();

        final DbPair dbLink = breakend.getDBLink();
        if(dbLink != null && dbLink.length() < 0)
        {
            cnLowSide -= dbLink.getOtherBreakend(breakend).jcn();
        }

        return cnLowSide;
    }

    public void findGermlineGeneDeletions(final List<SvCluster> clusters)
    {
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

                        for(GeneData geneData : mDisruptionGenes)
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
                                        breakend.getSV(), breakend.usesStart(), true, geneData, canonicalTrans,
                                        new int[] { 1, canonicalTrans.exons().size() + 1 }, UNKNOWN, UPSTREAM, 1.0);

                                mDisruptions.add(upDisruptionData);

                                SvDisruptionData downDisruptionData = new SvDisruptionData(
                                        nextBreakend.getSV(), nextBreakend.usesStart(), true, geneData, canonicalTrans,
                                        new int[] { 1, canonicalTrans.exons().size() + 1 }, UNKNOWN, DOWNSTREAM, 1.0);

                                mDisruptions.add(downDisruptionData);
                            }
                        }
                    }
                }
            }
        }
    }

    private static final String COHORT_WRITER_DISRUPTION = "Disruption";

    @Override
    public String fileType() { return COHORT_WRITER_DISRUPTION; }

    @Override
    public BufferedWriter createWriter(final String outputDir)
    {
        try
        {
            String outputFilename = outputDir + "LNX_DISRUPTIONS.csv";

            BufferedWriter writer = createBufferedWriter(outputFilename, false);

            writer.write("SampleId,Reportable,SvId,IsStart,Type,ClusterId,Chromosome,Position,Orientation");
            writer.write(",GeneId,GeneName,Strand,TransId,ExonUp,ExonDown,CodingType,RegionType");

            if(!mIsGermline)
            {
                writer.write(",UndisruptedCN,ExcludedReason,ExtraInfo");
            }
            else
            {
                writer.write(",ResolvedType,Filter,PonCount");
            }

            writer.newLine();
            return writer;
        }
        catch (final IOException e)
        {
            LNX_LOGGER.error("error writing disruptions: {}", e.toString());
            return null;
        }
    }

    public void writeMultiSampleData(final String sampleId, final List<SvVarData> svList)
    {
        List<String> outputLines = Lists.newArrayList();

        for(final SvDisruptionData disruptionData : mDisruptions)
        {
            final SvVarData var = disruptionData.Var;
            boolean isSvStart = disruptionData.IsStart;
            final GeneData gene = disruptionData.Gene;

            StringBuilder sb = new StringBuilder();

            sb.append(String.format("%s,%s,%d,%s,%s,%d,%s,%d,%d",
                    sampleId, disruptionData.Reportable, var.id(), isSvStart,
                    var != null ? var.type() : "", var != null ? var.getCluster().id() : -1,
                    disruptionData.Gene.Chromosome, var.position(isSvStart), var.orientation(isSvStart)));

            sb.append(String.format(",%s,%s,%d,%s,%d,%d,%s,%s",
                    gene.GeneId, gene.GeneName, gene.Strand, disruptionData.Transcript.TransName,
                    disruptionData.Exons[FS_UP], disruptionData.Exons[FS_DOWN], disruptionData.CodingType, disruptionData.RegionType));

            if(mIsGermline)
            {
                int ponCount = var.getSvData().filter().equals(PON_FILTER_PON) ? mGermlinePonCache.getPonCount(var) : 0;

                sb.append(String.format(",%s,%s,%d",
                        var.getCluster().getResolvedType(), var.getSvData().filter(), ponCount));
            }
            else
            {
                sb.append(String.format(",%.2f,,",disruptionData.UndisruptedCopyNumber));
            }

            outputLines.add(sb.toString());
        }

        if(!mIsGermline)
        {
            for(Map.Entry<BreakendTransData, String> entry : mRemovedDisruptions.entrySet())
            {
                final String exclusionInfo = entry.getValue();

                if(exclusionInfo.equals(NON_DISRUPT_REASON_SIMPLE_SV))
                    continue;

                final BreakendTransData transcript = entry.getKey();
                final BreakendGeneData gene = transcript.gene();
                final SvVarData var = svList.stream().filter(x -> x.id() == gene.id()).findFirst().orElse(null);

                StringBuilder sb = new StringBuilder();

                sb.append(String.format("%s,%s,%d,%s,%s,%d,%s,%d,%d",
                        sampleId, transcript.reportableDisruption(), gene.id(), gene.isStart(),
                        var != null ? var.type() : "", var != null ? var.getCluster().id() : -1,
                        gene.chromosome(), gene.position(), gene.orientation()));

                String exclusionReason = exclusionInfo;
                String extraInfo = "";

                String[] contextInfo = exclusionInfo.split(ITEM_DELIM);

                if(contextInfo.length == 2)
                {
                    exclusionReason = contextInfo[0];
                    extraInfo = contextInfo[1];
                }

                sb.append(String.format(",%s,%s,%d,%s,%d,%d,%s,%s,%.2f,%s,%s",
                        gene.StableId, gene.GeneName, gene.Strand, transcript.transName(),
                        transcript.ExonUpstream, transcript.ExonDownstream, transcript.codingType(),
                        transcript.regionType(), transcript.undisruptedCopyNumber(), exclusionReason, extraInfo));

                outputLines.add(sb.toString());
            }
        }

        mCohortDataWriter.write(this, outputLines);
    }
}
