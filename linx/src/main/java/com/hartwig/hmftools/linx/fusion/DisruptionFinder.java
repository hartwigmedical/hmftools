package com.hartwig.hmftools.linx.fusion;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.LinxOutput.ITEM_DELIM;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.types.ResolvedType.LINE;
import static com.hartwig.hmftools.linx.visualiser.file.VisGeneAnnotationType.DISRUPTION;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
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
import com.hartwig.hmftools.linx.germline.GermlineDisruptions;
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
    private final GermlineDisruptions mGermlineDisruptions;

    private final CohortDataWriter mCohortDataWriter;

    public static final int MAX_NON_DISRUPTED_CHAIN_LENGTH = 5000;

    public DisruptionFinder(
            final LinxConfig config, final EnsemblDataCache geneTransCache,
            final GermlinePonCache germlinePonCache, final CohortDataWriter cohortDataWriter, final VisSampleData visSampleData)
    {
        mGeneTransCache = geneTransCache;

        mIsGermline = config.IsGermline;
        mGermlineDisruptions = mIsGermline ? new GermlineDisruptions(config, geneTransCache, germlinePonCache) : null;

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

        Set<SvCluster> clusters = Sets.newHashSet();

        for(final SvVarData var : svList)
        {
            markTranscriptsDisruptive(var);

            // inferred SGLs are always non-disruptive
            if (var.isInferredSgl())
            {
                var.getGenesList(true).stream().forEach(x -> x.transcripts().stream().forEach(y -> y.setIsDisruptive(false)));
            }

            if(!var.getCluster().getChains().isEmpty())
                clusters.add(var.getCluster());
        }

        clusters.forEach(x -> x.getChains().forEach(y -> checkChainedBreakends(y)));
    }

    private static final String NON_DISRUPT_REASON_PARTIAL_DUP = "PartialDUP";
    private static final String NON_DISRUPT_REASON_SIMPLE_SV = "SimpleSV";
    private static final String NON_DISRUPT_REASON_LINE = "LINE";
    private static final String NON_DISRUPT_REASON_OTHER_NON_GENIC = "OtherNonGenic";
    private static final String NON_DISRUPT_REASON_SAME_INTRON = "SameIntron";
    private static final String NON_DISRUPT_REASON_TRAVERSED_SAME_INTRON = "TraversedSameIntron";

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

        if(var.isSimpleType())
        {
            markNonDisruptiveGeneTranscripts(genesStart, genesEnd, NON_DISRUPT_REASON_SIMPLE_SV);

            if(var.type() == DUP)
                markNonDisruptiveDup(var);
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

    private void markNonDisruptiveGeneTranscripts(
            final List<BreakendGeneData> genesStart, final List<BreakendGeneData> genesEnd, final String context)
    {
        // check for breakend transcripts from DELs and DUPs which start and end in the same intron
        for(final BreakendGeneData geneStart : genesStart)
        {
            final BreakendGeneData geneEnd = genesEnd.stream()
                    .filter(x -> x.StableId.equals(geneStart.StableId)).findFirst().orElse(null);

            if(geneEnd == null)
                continue;

            markNonDisruptiveTranscripts(geneStart.transcripts(), geneEnd.transcripts(), context);
        }
    }

    private boolean markNonDisruptiveTranscripts(final List<BreakendTransData> transList1, final List<BreakendTransData> transList2, final String context)
    {
        // look for matching transcripts which are both in the same non-exonic section
        boolean foundMatchingTrans = false;

        for (final BreakendTransData trans1 : transList1)
        {
            final BreakendTransData trans2 = transList2.stream()
                    .filter(x -> x.transId() == trans1.transId()).findFirst().orElse(null);

            if(trans2 == null)
                continue;

            if(markNonDisruptiveTranscript(trans1, trans2, context))
                foundMatchingTrans = true;
        }

        return foundMatchingTrans;
    }

    private boolean markNonDisruptiveTranscript(final BreakendTransData trans1, final BreakendTransData trans2, final String context)
    {
        if(trans1.ExonUpstream == trans2.ExonUpstream && !trans1.isExonic() && !trans2.isExonic())
        {
            markNonDisruptiveTranscript(trans1, context);
            markNonDisruptiveTranscript(trans2, context);
            return true;
        }

        return false;
    }

    private void markNonDisruptiveDup(final SvVarData var)
    {
        // DUPs are only disruptive if both ends are within the same transcript,
        // so need to test each end in turn for a matching transcript
        for(int se = SE_START; se <= SE_END; ++se)
        {
            List<BreakendGeneData> genesList = var.getGenesList(isStart(se));

            if(genesList.isEmpty())
                continue;

            List<BreakendGeneData> otherGenesList = var.getGenesList(!isStart(se));

            for(BreakendGeneData gene : genesList)
            {
                BreakendGeneData otherGene = otherGenesList.stream()
                        .filter(x -> x.StableId.equals(gene.StableId)).findFirst().orElse(null);

                if(otherGene == null)
                {
                    // other breakend isn't in any genic region
                    gene.transcripts().forEach(x -> markNonDisruptiveTranscript(x, NON_DISRUPT_REASON_PARTIAL_DUP));
                    continue;
                }
                else
                {
                    for(BreakendTransData trans : gene.transcripts())
                    {
                        BreakendTransData matchingTrans = otherGene.transcripts().stream()
                                .filter(x -> x.transId() == trans.transId()).findFirst().orElse(null);

                        if(matchingTrans == null)
                        {
                            markNonDisruptiveTranscript(trans, NON_DISRUPT_REASON_PARTIAL_DUP);
                        }
                        else
                        {
                            boolean withinTrans = trans.regionType() == EXONIC || trans.regionType() == INTRONIC;
                            boolean withinOtherTrans = matchingTrans.regionType() == EXONIC || matchingTrans.regionType() == INTRONIC;

                            if(!withinTrans || !withinOtherTrans)
                            {
                                markNonDisruptiveTranscript(trans, NON_DISRUPT_REASON_PARTIAL_DUP);
                                markNonDisruptiveTranscript(matchingTrans, NON_DISRUPT_REASON_PARTIAL_DUP);
                            }
                        }
                    }
                }

                if(gene.isUpstream())
                {
                    for(BreakendTransData upTrans : gene.transcripts())
                    {
                        final BreakendTransData downTrans = otherGene.transcripts().stream()
                                .filter(x -> x.transId() == upTrans.transId()).findFirst().orElse(null);

                        if(upTrans.ExonUpstream == 1 && downTrans.ExonDownstream <= 2 && !upTrans.isExonic())
                        {
                            markNonDisruptiveTranscript(upTrans, NON_DISRUPT_REASON_PARTIAL_DUP);
                            markNonDisruptiveTranscript(downTrans, NON_DISRUPT_REASON_PARTIAL_DUP);
                        }
                    }
                }
            }
        }
    }

    private void markNonDisruptiveTranscript(final BreakendTransData transcript, final String context)
    {
        if(!transcript.isDisruptive())
            return;

        transcript.setIsDisruptive(false);
        registerNonDisruptedTranscript(transcript, context);
    }

    private void registerNonDisruptedTranscript(final BreakendTransData transcript, final String context)
    {
        // keep track of non-disruptive breakends purely for logging purposes
        if(mIsGermline)
            return;

        if(mRemovedDisruptions.containsKey(transcript))
            return;

        if(!transcript.isCanonical() || !matchesDisruptionGene(transcript.gene()))
            return;

        LNX_LOGGER.debug("excluding gene({}) svId({}) reason({})", transcript.geneName(), transcript.gene().id(), context);

        mRemovedDisruptions.put(transcript, context);
    }

    // chaining checks
    private void checkChainedBreakends(final SvChain chain)
    {
        // unmark chain start & end if opposite orientation and within same intron
        // checkChainOpenBreakends(chain);

        // unmark a fully-intronic linked pair of breakends (INV or BND) if their other breakends are both non-genic
        checkOtherEndsNonGenic(chain);

        // unmark a fully-intronic linked pair of breakends if their other breakends are not already
        // disruptive DEL or DUP (implying other end is outside same intron)
        // checkSameIntronLinks(chain);

        // unmark breakend(s) if a) disruptive and b) follows a set of links which return to same intron in opposite orientation
        checkTraverseToSameIntron(chain);
    }

    private void checkChainOpenBreakends(final SvChain chain)
    {
        // the open ones are non-disruptive if they have opposite orientations and are both within the same intronic section
        SvBreakend chainStart = chain.getOpenBreakend(true);
        SvBreakend chainEnd = chain.getOpenBreakend(false);

        if(chainStart == null || chainEnd == null || chainStart.orientation() == chainEnd.orientation())
            return;

        List<BreakendGeneData> genesStart = chainStart.getGenesList();
        List<BreakendGeneData> genesEnd = chainEnd.getGenesList();

        if(genesStart.isEmpty() || genesEnd.isEmpty())
            return;

        for(BreakendGeneData geneStart : genesStart)
        {
            BreakendGeneData geneEnd = genesEnd.stream().filter(x -> x.GeneName.equals(geneStart.GeneName)).findFirst().orElse(null);

            if(geneEnd == null)
                continue;

            markNonDisruptiveTranscripts(geneStart.transcripts(), geneEnd.transcripts(), "ChainEnds");
        }
    }

    private void checkOtherEndsNonGenic(final SvChain chain)
    {
        for(final LinkedPair pair : chain.getLinkedPairs())
        {
            SvBreakend breakend1 = pair.firstBreakend();
            SvBreakend breakend2 = pair.secondBreakend();

            // have to ignore DELs since they can delete parts of the gene from outside it
            if(breakend1.getSV().type() == DEL || breakend2.getSV().type() == DEL)
                continue;

            List<BreakendGeneData> genes1 = breakend1.getGenesList();
            List<BreakendGeneData> genes2 = breakend2.getGenesList();

            if(genes1.isEmpty() || genes2.isEmpty())
                continue;

            List<BreakendTransData> transList1 = getDisruptedTranscripts(genes1);
            List<BreakendTransData> transList2 = getDisruptedTranscripts(genes2);

            if(transList1.isEmpty() || transList2.isEmpty())
                continue;

            SvBreakend otherBreakend1 = breakend1.getOtherBreakend();
            SvBreakend otherBreakend2 = breakend2.getOtherBreakend();

            boolean otherBreakend1NonGenic = otherBreakend1 != null ? otherBreakend1.getGenesList().isEmpty() : false;
            boolean otherBreakend2NonGenic = otherBreakend2 != null ? otherBreakend2.getGenesList().isEmpty() : false;

            if(otherBreakend1NonGenic && otherBreakend2NonGenic)
            {
                if(markNonDisruptiveTranscripts(transList1, transList2, NON_DISRUPT_REASON_OTHER_NON_GENIC))
                {
                    LNX_LOGGER.debug("pair({}) length({}) fully intronic)", pair, pair.length());
                }
            }
        }
    }

    private void checkSameIntronLinks(final SvChain chain)
    {
        // check each link in the chain, marking as non-disruptive any linked pair which is either fully contained with the same intron
        // or non-genic, and none can pass through another splice acceptor
        for(LinkedPair pair : chain.getLinkedPairs())
        {
            boolean traversesGene = pairTraversesGene(pair, 0, false);

            if(traversesGene)
                continue;

            SvBreakend breakend1 = pair.firstBreakend();
            SvBreakend breakend2 = pair.secondBreakend();

            List<BreakendGeneData> genes1 = breakend1.getGenesList();
            List<BreakendGeneData> genes2 = breakend2.getGenesList();

            if(genes1.isEmpty() || genes2.isEmpty())
                continue;

            List<BreakendTransData> transList1 = getDisruptedTranscripts(genes1);
            List<BreakendTransData> transList2 = getDisruptedTranscripts(genes2);

            if(transList1.isEmpty() || transList2.isEmpty())
                continue;

            // DELs and DUPs with disruptive breakends cannot then be marked non-disruptive just because they form an intronic TI
            if(breakend1.getSV().isSimpleType() || breakend2.getSV().isSimpleType())
                continue;

            for(BreakendGeneData gene1 : genes1)
            {
                BreakendGeneData gene2 = genes2.stream().filter(x -> x.GeneName.equals(gene1.GeneName)).findFirst().orElse(null);

                if(gene2 == null)
                    continue;

                markNonDisruptiveTranscripts(gene1.transcripts(), gene2.transcripts(), NON_DISRUPT_REASON_SAME_INTRON);
            }
        }
    }

    private void checkTraverseToSameIntron(final SvChain chain)
    {
        if(chain.getOpenBreakend(true) == null || chain.getOpenBreakend(false) == null)
            return;

        // walk through the chain from start to end, breakend by breakend
        int linkIndex = -1;
        // boolean useFirst = true;

        List<LinkedPair> links = chain.getLinkedPairs();

        while(true)
        {
            ++linkIndex;

            if(linkIndex >= links.size())
                break;

            LinkedPair pair = links.get(linkIndex);
            // SvBreakend breakend = useFirst ? pair.firstBreakend() : pair.secondBreakend();
            SvBreakend breakend = pair.firstBreakend().getOtherBreakend();

            List<BreakendGeneData> genes = breakend.getGenesList();
            if(genes.isEmpty())
                continue;

            List<BreakendTransData> transList = getDisruptedTranscripts(genes);
            if(transList.isEmpty())
                continue;

            // look ahead from this breakend to any which then return to this same intron with the opposite orientation
            int nextLinkIndex = linkIndex - 1; // so starts with the link from this breakend

            boolean foundSameIntronBreakend = false;
            final List<LinkedPair> traversedPairs = Lists.newArrayList();
            int chainLength = 0;

            while(true)
            {
                ++nextLinkIndex;

                if(nextLinkIndex >= links.size())
                    break;

                final LinkedPair nextPair = links.get(nextLinkIndex);

                // does this next link traverse another splice acceptor?
                boolean traversesGene = pairTraversesGene(nextPair, 0, false);

                if(traversesGene)
                    break;

                chainLength += nextPair.length();

                if(chainLength > MAX_NON_DISRUPTED_CHAIN_LENGTH)
                    break;

                // keep track of traversed pairs since they may be marked as non-disruptive too
                traversedPairs.add(nextPair);

                // does it return to the same intron and correct orientation for any transcripts
                final SvBreakend nextBreakend = nextPair.secondBreakend().getOtherBreakend();

                if(nextBreakend == null)
                    continue;

                if(nextBreakend.orientation() == breakend.orientation() || !nextBreakend.chromosome().equals(breakend.chromosome()))
                    continue;

                List<BreakendGeneData> otherGenes = nextBreakend.getGenesList();

                if(otherGenes.isEmpty() || otherGenes.stream().noneMatch(x -> genes.stream().anyMatch(y -> x.GeneName.equals(y.GeneName))))
                    continue;

                List<BreakendTransData> otherTransList = getDisruptedTranscripts(otherGenes);

                if(!otherTransList.isEmpty())
                {
                    if (markNonDisruptiveTranscripts(transList, otherTransList, NON_DISRUPT_REASON_SAME_INTRON))
                    {
                        foundSameIntronBreakend = true;

                        // mark any traversed pair's breakends non-disruptive regardless of their transcripts since in the context
                        // of this chain they aren't disruptive to other introns
                        for(LinkedPair traversedPair : traversedPairs)
                        {
                            SvBreakend breakend1 = traversedPair.firstBreakend();
                            SvBreakend breakend2 = traversedPair.secondBreakend();

                            List<BreakendTransData> traversedTrans1 = getDisruptedTranscripts(breakend1.getGenesList());
                            List<BreakendTransData> traversedTrans2 = getDisruptedTranscripts(breakend2.getGenesList());

                            if(!traversedTrans1.isEmpty() && !traversedTrans2.isEmpty())
                            {
                                markNonDisruptiveTranscripts(traversedTrans1, traversedTrans2, NON_DISRUPT_REASON_TRAVERSED_SAME_INTRON);
                            }
                        }

                        removeNonDisruptedTranscripts(transList);

                        LNX_LOGGER.debug("breakends({} & {}) chain({} links({} -> {}) return to same intron length({}) traversedPairs({})",
                                breakend, nextBreakend, chain.id(), linkIndex, nextLinkIndex, chainLength, traversedPairs.size());
                    }
                }

                break;
            }

            if(foundSameIntronBreakend)
                linkIndex = nextLinkIndex;
        }

        /*
        for(SvVarData var : chain.getSvList())
        {
            if(var.isSglBreakend())
                continue;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                SvBreakend breakend = var.getBreakend(se);

                List<BreakendGeneData> genes = var.getGenesList(breakend.usesStart());
                if(genes.isEmpty())
                    continue;

                List<BreakendTransData> transList = getDisruptedTranscripts(genes);
                if(transList.isEmpty())
                    continue;

                checkChainedTranscripts(breakend, chain, genes, transList);
            }
        }
        */
    }

    private void checkChainedTranscripts(
            final SvBreakend breakend, final SvChain chain, final List<BreakendGeneData> genes, final List<BreakendTransData> transList)
    {
        // check whether this breakend follows a chain which comes back to the same intron with correct orientation and
        // without traversing a splice acceptor. If so, any links it passes through are also marked as non-disruptive

        final List<LinkedPair> links = chain.getLinkedPairs();

        final List<LinkedPair> traversedPairs = Lists.newArrayList();

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

                // keep track of traversed pairs since they may be marked as non-disruptive too
                traversedPairs.add(nextPair);

                // does it return to the same intron and correct orientation for any transcripts
                final SvBreakend nextBreakend = traverseUp ?
                        nextPair.secondBreakend().getOtherBreakend() : nextPair.firstBreakend().getOtherBreakend();

                if(nextBreakend == null)
                    continue;

                if(nextBreakend.orientation() == breakend.orientation() || !nextBreakend.chromosome().equals(breakend.chromosome()))
                    continue;

                List<BreakendGeneData> otherGenes = nextBreakend.getGenesList();

                if(otherGenes.isEmpty() || otherGenes.stream().noneMatch(x -> genes.stream().anyMatch(y -> x.GeneName.equals(y.GeneName))))
                    continue;

                List<BreakendTransData> otherTransList = getDisruptedTranscripts(otherGenes);

                if(!otherTransList.isEmpty())
                {
                    if (markNonDisruptiveTranscripts(transList, otherTransList, NON_DISRUPT_REASON_SAME_INTRON))
                    {
                        // mark any traversed pair's breakends non-disruptive regardless of their transcripts since in the context
                        // of this chain they aren't disruptive to other introns
                        for(LinkedPair traversedPair : traversedPairs)
                        {
                            SvBreakend breakend1 = traversedPair.firstBreakend();
                            SvBreakend breakend2 = traversedPair.secondBreakend();

                            List<BreakendTransData> traversedTrans1 = getDisruptedTranscripts(breakend1.getGenesList());
                            List<BreakendTransData> traversedTrans2 = getDisruptedTranscripts(breakend2.getGenesList());

                            if(!traversedTrans1.isEmpty() && !traversedTrans2.isEmpty())
                            {
                                markNonDisruptiveTranscripts(traversedTrans1, traversedTrans2, NON_DISRUPT_REASON_TRAVERSED_SAME_INTRON);
                            }
                        }

                        removeNonDisruptedTranscripts(transList);

                        LNX_LOGGER.debug("breakends({} & {}) return to same intron, chain({}) links({}) length({}) traversedPairs({})",
                                breakend, nextBreakend, chain.id(), abs(index - startIndex), chainLength, traversedPairs.size());

                        if (transList.isEmpty())
                            return;
                    }
                }
            }
        }
    }

    public boolean pairTraversesGene(final LinkedPair pair, int fusionDirection, boolean isPrecodingUpstream)
    {
        // for this pair to not affect a fusion or be consider non-disruptive, the section it traverses cannot cross
        // any gene's splice acceptor with the same strand direction unless that is part of a fully traversed non-coding 5' exon

        // fusion-direction is the stream, or zero means ignore

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

    public void findReportableDisruptions(final List<SvVarData> svList, final List<SvCluster> clusters)
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
                                var,  gene.isStart(), gene.getGeneData(), transcript.TransData,
                                new int[] { transcript.ExonUpstream, transcript.ExonDownstream },
                                transcript.codingType(), transcript.regionType(), transcript.undisruptedCopyNumber());

                        disruptionData.setReportable(true);

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

        if(mIsGermline)
            mGermlineDisruptions.findGeneDeletions(clusters);
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

    // used by fusion logic only
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

            writer.write("SampleId,Reportable,SvId,IsStart,Type,Chromosome,Position,Orientation");
            writer.write(",GeneId,GeneName,Strand,TransId,ExonUp,ExonDown,CodingType,RegionType");
            writer.write(",ClusterId,ResolvedType,ClusterCount");

            if(!mIsGermline)
            {
                writer.write(",UndisruptedCN,ExcludedReason,ExtraInfo");
            }
            else
            {
                writer.write(",Filter,PonCount");
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

    public void writeCohortData(final String sampleId, final List<SvVarData> svList)
    {
        if(mIsGermline)
        {
            List<String> outputLines = mGermlineDisruptions.formCohortData(sampleId, mDisruptions);
            mCohortDataWriter.write(this, outputLines);
            return;
        }

        List<String> outputLines = Lists.newArrayList();

        for(final SvDisruptionData disruptionData : mDisruptions)
        {
            StringBuilder sb = new StringBuilder();

            sb.append(String.format("%s,%s", sampleId, disruptionData.asCsv()));

            sb.append(String.format(",%.2f,,",disruptionData.UndisruptedCopyNumber));

            outputLines.add(sb.toString());
        }

        for(Map.Entry<BreakendTransData, String> entry : mRemovedDisruptions.entrySet())
        {
            final String exclusionInfo = entry.getValue();

            if(exclusionInfo.equals(NON_DISRUPT_REASON_SIMPLE_SV))
                continue;

            final BreakendTransData transcript = entry.getKey();
            final BreakendGeneData gene = transcript.gene();
            final SvVarData var = svList.stream().filter(x -> x.id() == gene.id()).findFirst().orElse(null);

            if(var == null)
                continue;

            final SvCluster cluster = var.getCluster();

            StringBuilder sb = new StringBuilder();

            sb.append(String.format("%s,%s,%d,%s,%s,%s,%d,%d",
                    sampleId, transcript.reportableDisruption(), gene.id(), gene.isStart(),
                    var.type(), gene.chromosome(), gene.position(), gene.orientation()));

            String exclusionReason = exclusionInfo;
            String extraInfo = "";

            String[] contextInfo = exclusionInfo.split(ITEM_DELIM);

            if(contextInfo.length == 2)
            {
                exclusionReason = contextInfo[0];
                extraInfo = contextInfo[1];
            }

            sb.append(String.format(",%s,%s,%d,%s,%d,%d,%s,%s,%d,%s,%d",
                    gene.StableId, gene.GeneName, gene.Strand, transcript.transName(),
                    transcript.ExonUpstream, transcript.ExonDownstream, transcript.codingType(), transcript.regionType(),
                    cluster.id(), cluster.getResolvedType(), cluster.getSvCount()));

            sb.append(String.format(",%.2f,%s,%s", transcript.undisruptedCopyNumber(), exclusionReason, extraInfo));

            outputLines.add(sb.toString());
        }

        mCohortDataWriter.write(this, outputLines);
    }
}
