package com.hartwig.hmftools.linx.fusion;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.HOM_DEL_DISRUPTION;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.HOM_DUP_DISRUPTION;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.EXONIC;
import static com.hartwig.hmftools.common.gene.TranscriptRegionType.INTRONIC;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DEL;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.isStart;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.DUP;
import static com.hartwig.hmftools.linx.LinxConfig.LNX_LOGGER;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.formatJcn;
import static com.hartwig.hmftools.linx.annotators.PseudoGeneFinder.isPseudogeneDeletion;
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
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.linx.gene.BreakendGeneData;
import com.hartwig.hmftools.linx.gene.BreakendTransData;
import com.hartwig.hmftools.linx.CohortFileInterface;
import com.hartwig.hmftools.linx.LinxConfig;
import com.hartwig.hmftools.linx.CohortDataWriter;
import com.hartwig.hmftools.linx.chaining.SvChain;
import com.hartwig.hmftools.linx.germline.GermlineDisruptions;
import com.hartwig.hmftools.linx.types.DbPair;
import com.hartwig.hmftools.linx.types.LinkedPair;
import com.hartwig.hmftools.linx.types.SvBreakend;
import com.hartwig.hmftools.linx.types.SvCluster;
import com.hartwig.hmftools.linx.types.SvVarData;
import com.hartwig.hmftools.linx.visualiser.file.VisSampleData;

public class DisruptionFinder implements CohortFileInterface
{
    private final EnsemblDataCache mGeneTransCache;
    private final Map<String,List<String>> mDisruptionGeneTranscripts;
    private final List<DriverGene> mDriverGenes;
    private final VisSampleData mVisSampleData;

    private final List<SvDisruptionData> mDisruptions;
    private final Map<BreakendTransData,String> mRemovedDisruptions; // cached for diagnostic purposes

    private final boolean mIsGermline;
    private final GermlineDisruptions mGermlineDisruptions;

    private final CohortDataWriter mCohortDataWriter;

    public static final int MAX_NON_DISRUPTED_CHAIN_LENGTH = 5000;

    public DisruptionFinder(
            final LinxConfig config, final EnsemblDataCache geneTransCache,
            final CohortDataWriter cohortDataWriter, final VisSampleData visSampleData)
    {
        mGeneTransCache = geneTransCache;

        mIsGermline = config.IsGermline;
        mGermlineDisruptions = mIsGermline ? new GermlineDisruptions(config, geneTransCache) : null;
        mDriverGenes = config.DriverGenes;

        mDisruptionGeneTranscripts = getDisruptionGeneTranscripts(config.DriverGenes, !mIsGermline, geneTransCache);

        mVisSampleData = visSampleData;
        mCohortDataWriter = cohortDataWriter;

        mDisruptions = Lists.newArrayList();
        mRemovedDisruptions = Maps.newHashMap();
    }

    public static Map<String,List<String>> getDisruptionGeneTranscripts(
            final List<DriverGene> driverGenes, boolean onlyReportable, final EnsemblDataCache geneTransCache)
    {
        // builds a list of genes meeting reportable driver disruption criteria and includes any non-canonical transcript names
        Map<String,List<String>> geneTransMap = Maps.newHashMap();

        for(DriverGene driverGene : driverGenes)
        {
            if(onlyReportable && !driverGene.reportDisruption())
                continue;

            GeneData geneData = geneTransCache.getGeneDataByName(driverGene.gene());

            if(geneData != null)
            {
                geneTransMap.put(geneData.GeneId, driverGene.additionalReportedTranscripts());
            }
        }

        return geneTransMap;
    }

    public List<SvDisruptionData> getDisruptions() { return mDisruptions; }
    public GermlineDisruptions germlineDisruptions() { return mGermlineDisruptions; }

    private boolean matchesDisruptionGene(final String geneId)
    {
        return mDisruptionGeneTranscripts.containsKey(geneId);
    }

    public boolean matchesDisruptionTranscript(final String geneId, final TranscriptData transData)
    {
        if(!mDisruptionGeneTranscripts.containsKey(geneId))
            return false;

        if(transData.IsCanonical)
            return true;

        List<String> otherTrans = mDisruptionGeneTranscripts.get(geneId);
        return otherTrans.stream().anyMatch(x -> transData.TransName.equals(x));
    }

    public void addDisruptionGene(final String geneId)
    {
        mDisruptionGeneTranscripts.put(geneId, Lists.newArrayList());
    }

    public void markTranscriptsDisruptive(final List<SvVarData> svList)
    {
        mRemovedDisruptions.clear();

        Set<SvCluster> clusters = Sets.newHashSet();

        for(final SvVarData var : svList)
        {
            markTranscriptsDisruptive(var);

            // inferred SGLs are always non-disruptive
            if(var.isInferredSgl())
            {
                var.getGenesList(true).stream().forEach(x -> x.transcripts().stream().forEach(y -> y.setIsDisruptive(false)));
            }

            if(!var.getCluster().getChains().isEmpty())
                clusters.add(var.getCluster());
        }

        clusters.forEach(x -> x.getChains().forEach(y -> checkChainedBreakends(y)));

        // unmark a fully-intronic linked pair of breakends (INV or BND) if their other breakends are both non-genic
        clusters.forEach(x -> x.getChains().forEach(y -> checkOtherEndsNonGenic(y)));
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

        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(se == SE_END && var.isSglBreakend())
                continue;

            final List<BreakendGeneData> svGenes = se == SE_START ? genesStart : genesEnd;

            for(BreakendGeneData gene : svGenes)
            {
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

        for(final BreakendGeneData gene : genes)
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
                    .filter(x -> x.geneId().equals(geneStart.geneId())).findFirst().orElse(null);

            if(geneEnd == null)
                continue;

            markNonDisruptiveTranscripts(geneStart.transcripts(), geneEnd.transcripts(), context);
        }
    }

    private boolean markNonDisruptiveTranscripts(final List<BreakendTransData> transList1, final List<BreakendTransData> transList2, final String context)
    {
        // look for matching transcripts which are both in the same non-exonic section
        boolean foundMatchingTrans = false;

        for(final BreakendTransData trans1 : transList1)
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
        if(var.getCluster().getSvCount() > 1) // only apply these checks for simple (unclustered DUPs)
            return;

        // DUPs are only disruptive if both ends are within the same transcript or chained
        // so need to test each end in turn for a matching transcript
        for(int se = SE_START; se <= SE_END; ++se)
        {
            List<BreakendGeneData> genesList = var.getGenesList(isStart(se));

            if(genesList.isEmpty())
                continue;

            List<BreakendGeneData> otherBreakendGenes = var.getGenesList(!isStart(se));

            // ignore if the other breakend may be in a different gene (eg making a fusion)
            if(!otherBreakendGenes.isEmpty())
            {
                if(otherBreakendGenes.stream().anyMatch(x -> genesList.stream().anyMatch(y -> !y.geneName().equals(x.geneName()))))
                    continue;
            }

            for(BreakendGeneData gene : genesList)
            {
                if(otherBreakendGenes.isEmpty())
                {
                    // other breakend isn't in any genic region
                    gene.transcripts().forEach(x -> markNonDisruptiveTranscript(x, NON_DISRUPT_REASON_PARTIAL_DUP));
                    continue;
                }

                BreakendGeneData otherBreakendGeneData = otherBreakendGenes.stream()
                        .filter(x -> x.geneId().equals(gene.geneId())).findFirst().orElse(null);

                if(otherBreakendGeneData != null)
                {
                    // other breakend is in the same gene, so check each transcript for a matching transcript in the other breakend
                    for(BreakendTransData trans : gene.transcripts())
                    {
                        BreakendTransData matchingTrans = otherBreakendGeneData.transcripts().stream()
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

                    if(gene.isUpstream())
                    {
                        for(BreakendTransData upTrans : gene.transcripts())
                        {
                            final BreakendTransData downTrans = otherBreakendGeneData.transcripts().stream()
                                    .filter(x -> x.transId() == upTrans.transId()).findFirst().orElse(null);

                            if(downTrans != null && upTrans.ExonUpstream == 1 && downTrans.ExonDownstream <= 2 && !upTrans.isExonic())
                            {
                                markNonDisruptiveTranscript(upTrans, NON_DISRUPT_REASON_PARTIAL_DUP);
                                markNonDisruptiveTranscript(downTrans, NON_DISRUPT_REASON_PARTIAL_DUP);
                            }
                        }
                    }
                }
                else
                {
                    for(BreakendTransData trans : gene.transcripts())
                    {
                        markNonDisruptiveTranscript(trans, NON_DISRUPT_REASON_PARTIAL_DUP);

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

        if(!matchesDisruptionTranscript(transcript.gene().geneId(), transcript.TransData))
            return;

        LNX_LOGGER.debug("excluding gene({}) svId({}) reason({})", transcript.geneName(), transcript.gene().id(), context);

        mRemovedDisruptions.put(transcript, context);
    }

    // chaining checks
    private void checkChainedBreakends(final SvChain chain)
    {
        // unmark assembled links within a chain if the both the chain's ends are non-genic
        checkChainThroughIntrons(chain);

        // unmark breakend(s) if a) disruptive and b) follows a set of links which return to same intron in opposite orientation
        checkTraverseToSameIntron(chain);
    }

    private void checkChainThroughIntrons(final SvChain chain)
    {
        SvBreakend chainStart = chain.getOpenBreakend(true);
        SvBreakend chainEnd = chain.getOpenBreakend(false);

        if(chainStart == null || chainEnd == null)
            return;

        // have to ignore DELs since they can delete parts of the gene from outside it
        if(chainStart.getSV().type() == DEL || chainEnd.getSV().type() == DEL)
            return;

        List<BreakendGeneData> genesStart = chainStart.getGenesList();
        List<BreakendGeneData> genesEnd = chainEnd.getGenesList();

        // both chain ends must start in non-genic regions
        if(!genesStart.isEmpty() || !genesEnd.isEmpty())
        {
            if(genesStart.stream().anyMatch(x -> x.transcripts().stream().anyMatch(y -> y.isDisruptive())))
                return;

            if(genesEnd.stream().anyMatch(x -> x.transcripts().stream().anyMatch(y -> y.isDisruptive())))
                return;
        }

        for(final LinkedPair pair : chain.getLinkedPairs())
        {
            if(!pair.isAssembled())
                continue;

            SvBreakend breakend1 = pair.firstBreakend();
            SvBreakend breakend2 = pair.secondBreakend();

            if(breakend1.getSV().type() == DEL || breakend2.getSV().type() == DEL)
                continue;

            List<BreakendGeneData> genes1 = breakend1.getGenesList();
            List<BreakendGeneData> genes2 = breakend2.getGenesList();

            if(genes1.isEmpty() || genes2.isEmpty())
                continue;

            // neither breakend can be in a deletion bridge with a breakend in the same cluster and gene
            if(genes1.stream().anyMatch(x -> inDeletionBridgeWithinGene(chain, breakend1, x.geneName())))
                continue;
            else if(genes2.stream().anyMatch(x -> inDeletionBridgeWithinGene(chain, breakend2, x.geneName())))
                continue;

            for(BreakendGeneData gene1 : genes1)
            {
                BreakendGeneData gene2 = genes2.stream().filter(x -> x.geneId().equals(gene1.geneId())).findFirst().orElse(null);

                if(gene2 == null)
                    continue;

                for(BreakendTransData trans1 : gene1.transcripts())
                {
                    if(trans1.regionType() != INTRONIC)
                        continue;

                    // are the breakends within the same intro for any of the transcripts
                    BreakendTransData trans2 = gene2.transcripts().stream()
                            .filter(x -> x.transId() == trans1.transId())
                            .filter(x -> x.regionType() == INTRONIC)
                            .filter(x -> x.ExonDownstream == trans1.ExonDownstream)
                            .findFirst().orElse(null);

                    if(trans2 == null || (!trans1.isDisruptive() && !trans2.isDisruptive()))
                        continue;

                    if(markNonDisruptiveTranscript(trans1, trans2, NON_DISRUPT_REASON_SAME_INTRON))
                    {
                        LNX_LOGGER.debug("cluster({}) chain({}) pair({}) length({}) fully intronic)",
                                breakend1.getSV().getCluster().id(), chain.id(), pair, pair.baseLength());
                    }
                }
            }
        }
    }

    private void checkTraverseToSameIntron(final SvChain chain)
    {
        if(chain.getOpenBreakend(true) == null || chain.getOpenBreakend(false) == null)
            return;

        // walk through the chain from start to end, breakend by breakend
        int linkIndex = -1;

        List<LinkedPair> links = chain.getLinkedPairs();

        while(true)
        {
            ++linkIndex;

            if(linkIndex >= links.size())
                break;

            LinkedPair pair = links.get(linkIndex);
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

                chainLength += nextPair.baseLength();

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

                if(otherGenes.isEmpty() || otherGenes.stream().noneMatch(x -> genes.stream().anyMatch(y -> x.geneName().equals(y.geneName()))))
                    continue;

                List<BreakendTransData> otherTransList = getDisruptedTranscripts(otherGenes);

                if(!otherTransList.isEmpty())
                {
                    if(markNonDisruptiveTranscripts(transList, otherTransList, NON_DISRUPT_REASON_SAME_INTRON))
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
    }

    private boolean inDeletionBridgeWithinGene(final SvChain chain, final SvBreakend breakend, final String geneName)
    {
        // check if this breakend is in a deletion bridge with a disruptive breakend in the same gene but a different chain
        if(breakend.getDBLink() == null)
            return false;

        final SvBreakend otherBreakend = breakend.getDBLink().getOtherBreakend(breakend);

        if(otherBreakend.getCluster() != breakend.getCluster())
            return false;

        // breakends in the same chain can form DBs with each other and their TIs can be tested as usual
        if(otherBreakend.getCluster().findChain(otherBreakend.getSV()) == chain)
            return false;

        return otherBreakend.getGenesList().stream()
                .filter(x -> x.geneName().equals(geneName))
                .anyMatch(x -> x.transcripts().stream().anyMatch(y -> y.isDisruptive()));
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

            // neither breakend can be in a deletion bridge with a breakend in the same cluster and gene
            if(genes1.stream().anyMatch(x -> inDeletionBridgeWithinGene(chain, breakend1, x.geneName())))
                continue;
            else if(genes2.stream().anyMatch(x -> inDeletionBridgeWithinGene(chain, breakend2, x.geneName())))
                continue;

            SvBreakend otherBreakend1 = breakend1.getOtherBreakend();
            SvBreakend otherBreakend2 = breakend2.getOtherBreakend();

            boolean otherBreakend1NonGenic = otherBreakend1 != null ? otherBreakend1.getGenesList().isEmpty() : false;
            boolean otherBreakend2NonGenic = otherBreakend2 != null ? otherBreakend2.getGenesList().isEmpty() : false;

            if(otherBreakend1NonGenic && otherBreakend2NonGenic)
            {
                if(markNonDisruptiveTranscripts(transList1, transList2, NON_DISRUPT_REASON_OTHER_NON_GENIC))
                {
                    LNX_LOGGER.debug("pair({}) length({}) fully intronic)", pair, pair.positionDistance());
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
                    for(final ExonData exonData : transData.exons())
                    {
                        if(exonData.Rank == 1)
                            continue;

                        if((geneData.Strand == 1 && lowerPos <= exonData.Start && upperPos >= exonData.Start)
                        || (geneData.Strand == -1 && lowerPos <= exonData.End && upperPos >= exonData.End))
                        {
                            // allow an exon to be fully traversed if the upstream transcript is pre-coding
                            if(isPrecodingUpstream && lowerPos <= exonData.Start && upperPos >= exonData.End)
                            {
                                if(geneData.Strand == 1 && (transData.CodingStart == null || upperPos < transData.CodingStart))
                                    continue;
                                else if(geneData.Strand == -1 && (transData.CodingEnd == null || lowerPos > transData.CodingEnd))
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

        for(final SvVarData var : svList)
        {
            for(int be = SE_START; be <= SE_END; ++be)
            {
                if(be == SE_END && var.isSglBreakend())
                    continue;

                final List<BreakendGeneData> tsgGenesList = var.getGenesList(isStart(be)).stream()
                        .filter(x -> matchesDisruptionGene(x.geneId())).collect(Collectors.toList());

                if(tsgGenesList.isEmpty())
                    continue;

                for(BreakendGeneData gene : tsgGenesList)
                {
                    List<BreakendTransData> reportableDisruptions = gene.transcripts().stream()
                            .filter(BreakendTransData::isDisruptive)
                            .collect(Collectors.toList());

                    for(BreakendTransData transcript : reportableDisruptions)
                    {
                        if(!matchesDisruptionTranscript(gene.geneId(), transcript.TransData))
                            continue;

                        SvBreakend breakend = var.getBreakend(be);

                        if(!transcript.undisruptedCopyNumberSet())
                            transcript.setUndisruptedCopyNumber(getUndisruptedCopyNumber(breakend));

                        LNX_LOGGER.debug("var({}) breakend({}) gene({}) transcript({}) is disrupted, cnLowside({})",
                                var.id(), breakend, gene.geneName(), transcript.transName(), formatJcn(transcript.undisruptedCopyNumber()));

                        transcript.setReportableDisruption(true);

                        SvDisruptionData disruptionData = new SvDisruptionData(
                                var,  gene.isStart(), gene.GeneData, transcript.TransData,
                                new int[] { transcript.ExonUpstream, transcript.ExonDownstream },
                                transcript.codingType(), transcript.regionType(), transcript.undisruptedCopyNumber(),
                                breakend.copyNumber());

                        disruptionData.setReportable(true);

                        if(var.type() == DEL
                        && isPseudogeneDeletion(var, var.position(true), var.position(false), transcript.TransData))
                        {
                            disruptionData.markPseudogeneDeletion();
                        }

                        mDisruptions.add(disruptionData);

                        if(mVisSampleData != null)
                        {
                            mVisSampleData.addGeneExonData(
                                    var.getCluster().id(), gene.geneId(), gene.geneName(),
                                    "", 0, gene.chromosome(), DISRUPTION);
                        }
                    }
                }
            }
        }

        if(mIsGermline)
            mGermlineDisruptions.findGeneDeletions(clusters);
    }

    public void addReportableDisruptions(final List<DriverCatalog> driverCatalogs)
    {
        for(SvDisruptionData disruptionData : mDisruptions)
        {
            if(!disruptionData.reportable())
                continue;

            if(driverCatalogs.stream()
                    .filter(x -> x.gene().equals(disruptionData.Gene.GeneName))
                    .anyMatch(x -> x.driver() == HOM_DUP_DISRUPTION || x.driver() == HOM_DEL_DISRUPTION || x.driver() == DriverType.DISRUPTION))
            {
                continue;
            }

            DriverGene driverGene = mDriverGenes.stream().filter(x -> x.gene().equals(disruptionData.Gene.GeneName)).findFirst().orElse(null);

            DriverCatalog driverCatalog = ImmutableDriverCatalog.builder()
                    .driver(DriverType.DISRUPTION)
                    .category(driverGene != null ? driverGene.likelihoodType() : TSG)
                    .gene(disruptionData.Gene.GeneName)
                    .transcript(disruptionData.Transcript.TransName)
                    .isCanonical(disruptionData.Transcript.IsCanonical)
                    .chromosome(disruptionData.Gene.Chromosome)
                    .chromosomeBand(disruptionData.Gene.KaryotypeBand)
                    .likelihoodMethod(LikelihoodMethod.DISRUPTION)
                    .driverLikelihood(0)
                    .missense(0)
                    .nonsense(0)
                    .splice(0)
                    .inframe(0)
                    .frameshift(0)
                    .biallelic(disruptionData.UndisruptedCopyNumber < 0.5)
                    .minCopyNumber(disruptionData.UndisruptedCopyNumber)
                    .maxCopyNumber(disruptionData.MaxCopyNumber)
                    .build();

            driverCatalogs.add(driverCatalog);
        }
    }

    public void writeGermlineDisruptions(final String sampleId, final String outputDir)
    {
        mGermlineDisruptions.writeGermlineSVs(mDisruptions, sampleId, outputDir);
    }

    public static double getUndisruptedCopyNumber(final SvBreakend breakend)
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
                    .filter(x -> x.geneId().equals(transcript.GeneId))
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
                .filter(x -> x.geneId().equals(upTrans.gene().geneId())).findFirst().orElse(null);

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
                writer.write(GermlineDisruptions.csvHeader());
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
                    gene.geneId(), gene.geneName(), gene.strand(), transcript.transName(),
                    transcript.ExonUpstream, transcript.ExonDownstream, transcript.codingType(), transcript.regionType(),
                    cluster.id(), cluster.getResolvedType(), cluster.getSvCount()));

            sb.append(String.format(",%.2f,%s,%s", transcript.undisruptedCopyNumber(), exclusionReason, extraInfo));

            outputLines.add(sb.toString());
        }

        mCohortDataWriter.write(this, outputLines);
    }
}
