package com.hartwig.hmftools.linx.fusion_likelihood;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_LENGTHS;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getChromosomalArmLength;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.hasAnyPhaseMatch;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.haveOverlap;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.mapExonPhase;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.regionsPhaseMatched;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_NON_CODING;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_LONG_SAME_ARM;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_REMOTE;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_SHORT_INV;
import static com.hartwig.hmftools.linx.fusion_likelihood.LikelihoodCalc.calcOverlapBucketAreas;
import static com.hartwig.hmftools.linx.fusion_likelihood.LikelihoodCalc.setBucketLengthData;
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.checkAddCombinedGenePhaseRegion;
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.checkAddGenePhaseRegion;
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.mergePhaseRegions;
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.validateSimpleVsCombinedPhaseRegions;
import static com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection.nextTranscriptExons;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptExonData;
import com.hartwig.hmftools.linx.gene.SvGeneTranscriptCollection;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class CohortExpFusions
{
    // config
    private List<Long> mProximateBucketLengths;

    private final Map<String, List<GeneRangeData>> mChrForwardGeneDataMap;
    private final Map<String, List<GeneRangeData>> mChrReverseGeneDataMap;
    private Map<String, GeneRangeData> mGeneIdRangeDataMap;
    private Map<String,String> mGeneRangeDataMap; // populated from file

    private Map<String,Map<Integer,Long>> mDelGenePairCounts; // pair of gene-pairs to their overlap counts keyed by bucket index
    private Map<String,Map<Integer,Long>> mDupGenePairCounts;

    // global counts by type and buck length
    private List<Integer> mGlobalProximateCounts; // indexed as per the proximate lengths
    private int mGlobalShortInvCount;
    private int mGlobalLongDelDupInvCount;
    private long mArmLengthFactor;

    // cache the total counts (expressed as a length) for each unique phase combination per chromosome
    private Map<String, List<GenePhaseRegion>> mChrPhaseRegionsPosStrand;
    private Map<String, List<GenePhaseRegion>> mChrPhaseRegionsNegStrand;

    private boolean mLogVerbose;

    public static int PRE_GENE_3P_DISTANCE = 10000;
    public static int SHORT_INV_BUCKET = 100000;
    public static long MIN_BUCKET_LENGTH = 100;
    public static final double GENOME_BASE_COUNT = 3e9;
    public static final double MIN_FUSION_RATE = 1e-12;
    public static final String GENE_PAIR_DELIM = "_";

    private static final Logger LOGGER = LogManager.getLogger(CohortExpFusions.class);

    public CohortExpFusions()
    {
        mChrForwardGeneDataMap = Maps.newHashMap();
        mChrReverseGeneDataMap = Maps.newHashMap();
        mGeneIdRangeDataMap = null;
        mDelGenePairCounts = Maps.newHashMap();
        mDupGenePairCounts = Maps.newHashMap();
        mChrPhaseRegionsPosStrand = Maps.newHashMap();
        mChrPhaseRegionsNegStrand = Maps.newHashMap();
        mGlobalProximateCounts = Lists.newArrayList();
        mGlobalShortInvCount = 0;
        mGlobalLongDelDupInvCount = 0;
        mArmLengthFactor = 0;
    }

    public final Map<String, List<GeneRangeData>> getChrGeneRangeDataMap() { return mChrForwardGeneDataMap; }
    public Map<String,Map<Integer,Long>> getDelGenePairCounts() { return mDelGenePairCounts; }
    public Map<String,Map<Integer,Long>> getDupGenePairCounts() { return mDupGenePairCounts; }

    public long getArmLengthFactor() { return mArmLengthFactor; }

    public void initialise(List<Long> proximateBucketLengths, int shortInvBucketLength, boolean logVerbose)
    {
        mProximateBucketLengths = proximateBucketLengths;

        mProximateBucketLengths.stream().forEach(x -> mGlobalProximateCounts.add(0));

        mLogVerbose = logVerbose;

        // sum up all arm lengths to adjusted same-arm fusion rates
        long maxBucketLength = max(getMaxBucketLength(), 4000000);

        for(Map.Entry<String,Integer> entry : CHROMOSOME_LENGTHS.entrySet())
        {
            final String chromosome = entry.getKey();

            long armLength = getChromosomalArmLength(chromosome, CHROMOSOME_ARM_P);

            if(armLength > maxBucketLength)
                mArmLengthFactor += pow(armLength - maxBucketLength, 2);

            armLength = getChromosomalArmLength(chromosome, CHROMOSOME_ARM_Q);

            if(armLength > maxBucketLength)
                mArmLengthFactor += pow(armLength - maxBucketLength, 2);
        }
    }

    public void initialiseGeneIdRangeDataMap()
    {
        mGeneIdRangeDataMap = Maps.newHashMap();
    }

    public GeneRangeData findGeneRangeData(final String geneId)
    {
        return mGeneIdRangeDataMap != null ? mGeneIdRangeDataMap.get(geneId) : null;
    }

    public void generateGenePhasingCounts(final SvGeneTranscriptCollection geneTransCache,
            List<String> restrictedChromosomes, List<String> restrictedGeneIds)
    {
        // for each gene, walk through all its transcripts and count up the number of bases in each phasing region (eg 5'UTR, 0-2)
        // also convert these into combined phasing regions (eg where say phase 0 and 1 overlap) for using in downstream analysis
        // finally gather up the total bases by phasing type across each chromosome
        final Map<String, List<EnsemblGeneData>> chrGeneDataMap = geneTransCache.getChrGeneDataMap();
        final Map<String, List<EnsemblGeneData>> chrReverseGeneDataMap = geneTransCache.getChrReverseGeneDataMap();

        for(Map.Entry<String, List<EnsemblGeneData>> entry : chrGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            if(!restrictedChromosomes.isEmpty() && !restrictedChromosomes.contains(chromosome))
                continue;

            LOGGER.info("generating phase counts for chromosome({})", chromosome);

            final List<EnsemblGeneData> geneDataList = entry.getValue();
            final List<EnsemblGeneData> reverseGeneDataList = chrReverseGeneDataMap.get(chromosome);

            List<GeneRangeData> geneList = Lists.newArrayList();
            List<GeneRangeData> geneEndFirstList = Lists.newArrayList();
            List<GenePhaseRegion> chrPhaseRegionsPosStrand = Lists.newArrayList();
            List<GenePhaseRegion> chrPhaseRegionsNegStrand = Lists.newArrayList();
            long posStrandTotal = 0;
            long negStrandTotal = 0;

            for(final EnsemblGeneData geneData :entry.getValue())
            {
                GeneRangeData geneRangeData = new GeneRangeData(geneData);

                geneList.add(geneRangeData);

                int index = 0;
                for(; index < geneEndFirstList.size(); ++index)
                {
                    final GeneRangeData rgd = geneEndFirstList.get(index);

                    if(geneData.GeneEnd < rgd.GeneData.GeneEnd)
                        break;
                }

                geneEndFirstList.add(index, geneRangeData);

                if(!restrictedGeneIds.isEmpty() && !restrictedGeneIds.contains(geneData.GeneId))
                    continue;

                // load from Ensembl transcript and exon data
                final List<TranscriptExonData> transExonDataList = geneTransCache.getTransExonData(geneData.GeneId);

                if (transExonDataList == null)
                    continue;

                long precedingGeneSAPos = geneTransCache.findPrecedingGeneSpliceAcceptorPosition(
                        geneData, geneData.Strand == 1 ? reverseGeneDataList : geneDataList);

                if(geneData.GeneId.equals("ENSG00000260411"))
                {
                    LOGGER.info("spec gene({})", geneData.GeneId);
                }

                generateGenePhaseRegions(geneRangeData, transExonDataList, precedingGeneSAPos);

                List<GenePhaseRegion> combinedPhaseRegions = geneRangeData.getPhaseRegions();

                /*
                // convert to a non-overlapping combined set of phase regions
                combinedPhaseRegions = Lists.newArrayList();

                for (GenePhaseRegion region : geneRangeData.getPhaseRegions())
                {
                    if (region.length() <= 0)
                    {
                        LOGGER.warn("gene({}) invalid length({})", region.GeneId, region.length());
                        continue;
                    }

                    checkAddCombinedGenePhaseRegion(region, combinedPhaseRegions);
                }
                */

                geneRangeData.setCombinedPhaseRegions(combinedPhaseRegions);
                mergePhaseRegions(combinedPhaseRegions);

                validateSimpleVsCombinedPhaseRegions(geneData.GeneId, geneRangeData.getPhaseRegions(), combinedPhaseRegions);

                if(mLogVerbose)
                {
                    for(GenePhaseRegion region : geneRangeData.getPhaseRegions())
                    {
                        LOGGER.debug("gene({}) range({} -> {}) phase({}) preGene({})",
                                region.GeneId, region.start(), region.end(), region.Phase, region.isAnyPreGene());
                    }
                }

                for(GenePhaseRegion region : combinedPhaseRegions)
                {
                    // validity check
                    if(region.length() <= 0)
                    {
                        LOGGER.error("invalid region: gene({}) range({} -> {}) phase({})",
                                region.GeneId, region.start(), region.end(), region.getCombinedPhase());
                        continue;
                    }

                    // add to the chromosome's summary phasing regions by strand, without concern for overlaps
                    // across genes on the same strand
                    if(geneData.Strand == 1)
                    {
                        posStrandTotal += region.length();
                        addPhaseRegion(chrPhaseRegionsPosStrand, region);
                    }
                    else
                    {
                        negStrandTotal += region.length();
                        addPhaseRegion(chrPhaseRegionsNegStrand, region);
                    }
                }
            }

            mChrPhaseRegionsPosStrand.put(chromosome, chrPhaseRegionsPosStrand);
            mChrPhaseRegionsNegStrand.put(chromosome, chrPhaseRegionsNegStrand);

            LOGGER.debug("chr({}) posStrandTotal({}) negStrandTotal({})", chromosome, posStrandTotal, negStrandTotal);

            mChrForwardGeneDataMap.put(chromosome, geneList);
            mChrReverseGeneDataMap.put(chromosome, geneEndFirstList);
        }
    }

    private void addPhaseRegion(List<GenePhaseRegion> regions, final GenePhaseRegion newRegion)
    {
        // add if new, else expand existing
        for(GenePhaseRegion region : regions)
        {
            if(region.getCombinedPhase() == newRegion.getCombinedPhase())
            {
                region.setEnd(region.end() + newRegion.length());
                return;
            }
        }

        regions.add(new GenePhaseRegion("", 0, newRegion.length(),
                newRegion.getPhaseArray(), newRegion.getPreGenePhaseStatus()));
    }

    public void generateGenePhaseRegions(
            GeneRangeData geneRangeData, final List<TranscriptExonData> transExonDataList, long precedingGeneSAPos)
    {
        // 5-prime rules - must be post-promotor
        // 3-prime rules: must be coding and > 1 exon, needs to find first splice acceptor and then uses its phasing
        // need to mark coding and non-coding regions
        List<GenePhaseRegion> phaseRegions = Lists.newArrayList();

        final EnsemblGeneData geneData = geneRangeData.GeneData;

        List<GenePhaseRegion> matchedTransRegions = Lists.newArrayList();

        int teIndex = 0;
        List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

        while (!transcriptExons.isEmpty())
        {
            List<GenePhaseRegion> transcriptRegions = Lists.newArrayList();

            for (int i = 0; i < transcriptExons.size() - 1; ++i)
            {
                // skip single-exon transcripts since without an intronic section their fusion likelihood is negligible
                if(transcriptExons.size() == 1)
                    break;

                TranscriptExonData exonData = transcriptExons.get(i);

                if (exonData.CodingStart == null)
                {
                    // mark the whole transcript as a single UTR
                    GenePhaseRegion phaseRegion = new GenePhaseRegion(
                            geneData.GeneId, exonData.TransStart, exonData.TransEnd, PHASE_NON_CODING);

                    checkAddCombinedGenePhaseRegion(phaseRegion, phaseRegions);

                    GenePhaseRegion transRegion = GenePhaseRegion.from(phaseRegion);
                    transcriptRegions.add(transRegion);
                    break;
                }

                TranscriptExonData nextExonData = transcriptExons.get(i+1);

                if(geneData.Strand == 1 && nextExonData.ExonStart > exonData.CodingEnd)
                    break;
                else if(geneData.Strand == -1 && exonData.ExonEnd < exonData.CodingStart)
                    continue;

                // add an upstream gene region (only applicable for downstream fusion genes)
                // with a phase of -1 unless coding starts in the first exon or on the first base of the second exon
                if(precedingGeneSAPos > 0 && geneData.Strand == 1 && i == 0)
                {
                    int regionPhase = (exonData.CodingStart < exonData.ExonEnd) ? nextExonData.ExonPhase : -1;

                    long preDistance = max(exonData.TransStart - precedingGeneSAPos, 0);
                    long upstreamDistance = min(preDistance, PRE_GENE_3P_DISTANCE);
                    long regionStart = exonData.TransStart - upstreamDistance;
                    long regionEnd = exonData.ExonStart - 1;

                    if(regionStart < regionEnd)
                    {
                        GenePhaseRegion phaseRegion =
                                new GenePhaseRegion(geneData.GeneId, regionStart, regionEnd, mapExonPhase(regionPhase));
                        phaseRegion.setPreGene(true, phaseRegion.Phase);


                        checkAddCombinedGenePhaseRegion(phaseRegion, phaseRegions);
                    }
                }
                else if(precedingGeneSAPos > 0 && geneData.Strand == -1 && i == transcriptExons.size() - 2)
                {
                    int regionPhase = (exonData.CodingEnd > nextExonData.ExonStart) ? exonData.ExonPhase : -1;

                    long regionStart = nextExonData.ExonEnd + 1;

                    long preDistance = max(precedingGeneSAPos - exonData.TransEnd, 0);
                    long upstreamDistance = min(preDistance, PRE_GENE_3P_DISTANCE);
                    long regionEnd = nextExonData.TransEnd + upstreamDistance;

                    if(regionStart < regionEnd)
                    {
                        GenePhaseRegion phaseRegion =
                                new GenePhaseRegion(geneData.GeneId, regionStart, regionEnd, mapExonPhase(regionPhase));
                        phaseRegion.setPreGene(true, phaseRegion.Phase);

                        checkAddCombinedGenePhaseRegion(phaseRegion, phaseRegions);

                        GenePhaseRegion transRegion = GenePhaseRegion.from(phaseRegion);
                        transcriptRegions.add(transRegion);
                    }
                }

                // turn the intronic section into a phase region and fold the exon in with same phasing for simplicity
                // step back a single base so the consecutive regions do not overlap

                GenePhaseRegion phaseRegion = null;
                if(geneData.Strand == 1)
                {
                    phaseRegion = new GenePhaseRegion(
                            geneData.GeneId, exonData.ExonStart, nextExonData.ExonStart - 1, mapExonPhase(exonData.ExonPhaseEnd));
                }
                else
                {
                    phaseRegion = new GenePhaseRegion(
                            geneData.GeneId, exonData.ExonEnd + 1, nextExonData.ExonEnd, mapExonPhase(nextExonData.ExonPhaseEnd));
                }

                checkAddCombinedGenePhaseRegion(phaseRegion, phaseRegions);

                GenePhaseRegion transRegion = GenePhaseRegion.from(phaseRegion);
                transcriptRegions.add(transRegion);
            }

            generateSameGeneCounts(geneRangeData, transcriptRegions, matchedTransRegions);

            teIndex += transcriptExons.size();
            transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
        }

        geneRangeData.setPhaseRegions(phaseRegions);

    }

    public void generateSameGeneCounts(GeneRangeData geneRangeData, final List<GenePhaseRegion> transcriptRegions,
            List<GenePhaseRegion> matchedTransRegions)
    {
        // look for same-gene fusions - first removing or shrinking any regions from this gene already used

        // first reduce or remove any regions which overlap with those already causing a same-gene fusion
        if(!matchedTransRegions.isEmpty())
        {
            int index = 0;
            while (index < transcriptRegions.size())
            {
                GenePhaseRegion newRegion = transcriptRegions.get(index);

                boolean fullyCovered = false;
                for (final GenePhaseRegion region : matchedTransRegions)
                {
                    if (!haveOverlap(region, newRegion, false))
                        continue;

                    if (region.start() <= newRegion.start() && region.end() >= newRegion.end())
                    {
                        // fully covered so remove
                        newRegion.setEnd(newRegion.start());
                        break;
                    }
                    else if (newRegion.start() <= region.start() && newRegion.end() >= region.end())
                    {
                        // split and add a new region after the existing
                        GenePhaseRegion extraRegion = GenePhaseRegion.from(newRegion);
                        extraRegion.setStart(region.end() + 1);
                        transcriptRegions.add(extraRegion);

                        newRegion.setEnd(region.start() - 1);
                    }
                    else if (newRegion.start() < region.start())
                    {
                        newRegion.setEnd(region.start() - 1);
                    }
                    else if (newRegion.end() > region.end())
                    {
                        newRegion.setStart(region.end() + 1);
                    }
                }

                if (newRegion.length() <= 1)
                    transcriptRegions.remove(index);
                else
                    ++index;
            }
        }

        // now look for candidate fusions within these transcript's remaining regions
        for (int i = 0; i < transcriptRegions.size(); ++i)
        {
            GenePhaseRegion region1 = transcriptRegions.get(i);

            for (int j = i+1; j < transcriptRegions.size(); ++j)
            {
                GenePhaseRegion region2 = transcriptRegions.get(j);
                if(testProximatePhaseRegions(geneRangeData, geneRangeData, region1, region2))
                {
                    if(!matchedTransRegions.contains(region1))
                        matchedTransRegions.add(region1);

                    if(!matchedTransRegions.contains(region2))
                        matchedTransRegions.add(region2);
                }
            }
        }
    }

    public void generateNonProximateCounts()
    {
        // for each arm and each gene in that arm sum up the overlapping counts against all genes beyond the specified
        // DEL and DUP max bucket length, and then all overlapping counts on remote arms
        long proximateLimit = getMaxBucketLength();

        for (Map.Entry<String, List<GeneRangeData>> entry : mChrForwardGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            LOGGER.info("calculating local and remote overlap counts for chr({})", chromosome);

            List<GeneRangeData> geneList = entry.getValue();

            for(int i = 0; i < geneList.size(); ++i)
            {
                GeneRangeData gene1 = geneList.get(i);

                /*
                if(gene1.GeneData.GeneId.equals("ENSG00000189283"))
                {
                    LOGGER.info("specific gene: {}", gene1.GeneData.GeneId);
                }
                */

                // determine the remote overlap count by comparing this gene's regions to all matching remote ones
                for(int j = 0; j <= 1; ++j)
                {
                    final Map<String, List<GenePhaseRegion>> chrPhaseRegions =
                            (j == 0) ? mChrPhaseRegionsPosStrand : mChrPhaseRegionsNegStrand;

                    for (Map.Entry<String, List<GenePhaseRegion>> chrEntry : chrPhaseRegions.entrySet())
                    {
                        if (chrEntry.getKey().equals(chromosome)) // skip the same chromosome
                            continue;

                        final List<GenePhaseRegion> remoteRegions = chrEntry.getValue();

                        for (GenePhaseRegion region : gene1.getCombinedPhaseRegions())
                        {
                            for (GenePhaseRegion remoteRegion : remoteRegions)
                            {
                                // test gene as an downstream vs all remote upstream phasing regions
                                if (hasAnyPhaseMatch(remoteRegion, region, false) || regionsPhaseMatched(remoteRegion, region))
                                {
                                    gene1.addBaseOverlapCountDownstream(NON_PROX_TYPE_REMOTE, region.length() * remoteRegion.length());
                                }

                                // then as an upstream partner
                                if (hasAnyPhaseMatch(region, remoteRegion, false) || regionsPhaseMatched(region, remoteRegion))
                                {
                                    gene1.addBaseOverlapCountUpstream(NON_PROX_TYPE_REMOTE, region.length() * remoteRegion.length());
                                }
                            }
                        }
                    }
                }

                // and now the local ones outside the DEL and DUP proximity lengths, ignoring any strand checking
                // this is considering fusions between DUPs, DUPs or INVs
                for(int j = i; j < geneList.size(); ++j)
                {
                    GeneRangeData gene2 = geneList.get(j);

                    if(!gene1.Arm.equals(gene2.Arm))
                        break;

                    // find closest and furthest distance between these 2 genes
                    long minDistance = min(abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneStart),
                            abs(gene1.GeneData.GeneEnd - gene2.GeneData.GeneEnd));

                    minDistance = min(minDistance, abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneEnd));
                    minDistance = min(minDistance, abs(gene1.GeneData.GeneEnd - gene2.GeneData.GeneStart));

                    long maxDistance = max(abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneStart),
                            abs(gene1.GeneData.GeneEnd - gene2.GeneData.GeneEnd));

                    maxDistance = max(maxDistance, abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneEnd));
                    maxDistance = max(maxDistance, abs(gene1.GeneData.GeneEnd - gene2.GeneData.GeneStart));

                    int type = -1;

                    // assign the proximity to a category to build up global stats
                    if (j > i && minDistance > proximateLimit)
                    {
                        ++mGlobalLongDelDupInvCount;
                        type = NON_PROX_TYPE_LONG_SAME_ARM;
                    }
                    else if(gene1.GeneData.Strand != gene2.GeneData.Strand && maxDistance < SHORT_INV_BUCKET)
                    {
                        ++mGlobalShortInvCount;
                        type = NON_PROX_TYPE_SHORT_INV;
                    }
                    else
                    {
                        // assign proximate genes to correct bucket length, same-genes will be counted here as well
                        if(gene1.GeneData.Strand == gene2.GeneData.Strand)
                        {
                            for (int b = 0; b < mProximateBucketLengths.size() - 1; ++b)
                            {
                                long minLength = mProximateBucketLengths.get(b);
                                long maxLength = mProximateBucketLengths.get(b + 1);

                                if (minDistance > maxLength || maxDistance < minLength)
                                    continue;

                                mGlobalProximateCounts.set(b, mGlobalProximateCounts.get(b) + 1);
                            }
                        }

                        continue;
                    }

                    // calculate phasing overlap areas
                    for (GenePhaseRegion region1 : gene1.getCombinedPhaseRegions())
                    {
                        // the downstream gene of the potential fusion cannot be non-coding
                        for (GenePhaseRegion region2 : gene2.getCombinedPhaseRegions())
                        {
                            long regionOverlap = region1.length() * region2.length();

                            if(hasAnyPhaseMatch(region1, region2, false))
                            {
                                gene1.addBaseOverlapCountUpstream(type, regionOverlap);
                                gene1.addBaseOverlapCountDownstream(type, regionOverlap);

                                gene2.addBaseOverlapCountUpstream(type, regionOverlap);
                                gene2.addBaseOverlapCountDownstream(type, regionOverlap);
                            }
                            else
                            {
                                // check each direction in turn
                                if (regionsPhaseMatched(region1, region2))
                                {
                                    gene1.addBaseOverlapCountUpstream(type, regionOverlap);
                                    gene2.addBaseOverlapCountDownstream(type, regionOverlap);
                                }

                                if (regionsPhaseMatched(region2, region1))
                                {
                                    gene2.addBaseOverlapCountUpstream(type, regionOverlap);
                                    gene1.addBaseOverlapCountDownstream(type, regionOverlap);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public void generateProximateFusionCounts()
    {
        if(mProximateBucketLengths.isEmpty())
            return;

        // walk through each chromosome, and test distances between each gene and the next ones for potential
        // fusions from DELs and DUPs within the configured bucket distances

        LOGGER.info("finding proximate fusion candidates on forward strand");

        for (Map.Entry<String, List<GeneRangeData>> entry : mChrForwardGeneDataMap.entrySet())
        {
            LOGGER.info("proximate forward-strand fusion candidates from chromosome({})", entry.getKey());
            findProximateFusions(entry.getValue(), 1);
        }

        LOGGER.info("finding proximate fusion candidates on reverse strand");

        for (Map.Entry<String, List<GeneRangeData>> entry : mChrReverseGeneDataMap.entrySet())
        {
            LOGGER.info("proximate reverse-strand fusion candidates from chromosome({})", entry.getKey());
            findProximateFusions(entry.getValue(), -1);
        }
    }

    public long getMaxBucketLength()
    {
        return !mProximateBucketLengths.isEmpty() ? mProximateBucketLengths.get(mProximateBucketLengths.size() - 1) : MIN_BUCKET_LENGTH;
    }

    private void findProximateFusions(List<GeneRangeData> geneRangeList, int strandMatch)
    {
        long rangeLimit = getMaxBucketLength();
        long minBucketLength = mProximateBucketLengths.get(0);

        for(int lowerIndex = 0; lowerIndex < geneRangeList.size(); ++lowerIndex)
        {
            GeneRangeData lowerGene = geneRangeList.get(lowerIndex);

            if(lowerGene.GeneData.Strand != strandMatch)
                continue;

            if(lowerGene.GeneData.GeneId.equals("ENSG00000260411"))
            {
                LOGGER.info("spec gene({})", lowerGene.GeneData.GeneId);
            }

            // don't allow same-gene fusions (they are handled within a transcript), so start the index at the next gene
            for(int upperIndex = lowerIndex+1; upperIndex < geneRangeList.size(); ++upperIndex)
            {
                GeneRangeData upperGene = geneRangeList.get(upperIndex);

                if(upperGene.GeneData.Strand != strandMatch)
                    continue;

                // exit if the gene is now too far away from the lower gene
                if(upperGene.GeneData.GeneStart - lowerGene.GeneData.GeneEnd > rangeLimit)
                    break;

                // find all matching phasing regions and account for any duplicated overlapping regions from different phasing matches
                // for all those found, assign them to a length bucket and find the number of bases that fall into the length bucket region
                for (GenePhaseRegion lowerRegion : lowerGene.getCombinedPhaseRegions())
                {
                    for (GenePhaseRegion upperRegion : upperGene.getCombinedPhaseRegions())
                    {
                        if (!haveOverlap(lowerRegion, upperRegion, false))
                        {
                            testProximatePhaseRegions(lowerGene, upperGene, lowerRegion, upperRegion);
                        }
                        else
                        {
                            // first process the regions before and after the overlapping region
                            GenePhaseRegion lowerSplit = new GenePhaseRegion(lowerRegion.GeneId, lowerRegion.start(), upperRegion.start() - 1,
                                    lowerRegion.getPhaseArray(), lowerRegion.getPreGenePhaseStatus());

                            testProximatePhaseRegions(lowerGene, upperGene, lowerSplit, upperRegion);

                            if(lowerRegion.end() > upperRegion.end())
                            {
                                GenePhaseRegion upperSplit = new GenePhaseRegion(lowerRegion.GeneId, upperRegion.end() + 1, lowerRegion.end(),
                                        lowerRegion.getPhaseArray(), lowerRegion.getPreGenePhaseStatus());

                                testProximatePhaseRegions(upperGene, lowerGene, upperRegion, upperSplit);
                            }
                            else
                            {
                                GenePhaseRegion upperSplit = new GenePhaseRegion(upperRegion.GeneId, lowerRegion.end() + 1, upperRegion.end(),
                                        upperRegion.getPhaseArray(), upperRegion.getPreGenePhaseStatus());

                                testProximatePhaseRegions(lowerGene, upperGene, lowerRegion, upperSplit);
                            }

                            // next split the overlap into smaller chunks to proces these one by one
                            long overlapStart = upperRegion.start();
                            long overlapEnd = min(lowerRegion.end(), upperRegion.end());

                            double segment = (overlapEnd - overlapStart) / 5.0;
                            long segmentSize = max((long)floor(segment), minBucketLength);

                            // break the regions up into non-overlapping pieces and process them one by one
                            List<GenePhaseRegion> lowerSplits = Lists.newArrayList();
                            List<GenePhaseRegion> upperSplits = Lists.newArrayList();

                            for(int i = 0; i < 5; ++i)
                            {
                                long regionStart = overlapStart + (i * segmentSize);
                                long regionEnd = regionStart + segmentSize - 1;

                                if(regionEnd > overlapEnd)
                                    break;

                                lowerSplits.add(new GenePhaseRegion(lowerRegion.GeneId, regionStart, regionEnd,
                                        lowerRegion.getPhaseArray(), lowerRegion.getPreGenePhaseStatus()));

                                upperSplits.add(new GenePhaseRegion(upperRegion.GeneId, regionStart, regionEnd,
                                        upperRegion.getPhaseArray(), upperRegion.getPreGenePhaseStatus()));
                            }

                            for(int i = 0; i < 5; ++i)
                            {
                                GenePhaseRegion lowerSegment = lowerSplits.get(i);

                                for(int j = i + 1; j < 5; ++j)
                                {
                                    GenePhaseRegion upperSegment = upperSplits.get(j);
                                    testProximatePhaseRegions(lowerGene, upperGene, lowerSegment, upperSegment);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public boolean testProximatePhaseRegions(GeneRangeData gene1, GeneRangeData gene2, GenePhaseRegion region1, GenePhaseRegion region2)
    {
        // ignore overlapping regions for now since it's not clear whether a DUP or DEL would be required
        if (haveOverlap(region1, region2, false))
            return false;

        // test for DEL and DUP fusion independently
        boolean isForwardStrand = (gene1.GeneData.Strand == 1);
        boolean foundMatch = false;

        for(int i = 0; i <= 1; ++i)
        {
            boolean isDel = (i == 0);

            boolean region1IsLower = region1.start() < region2.start();
            GeneRangeData lowerGene = region1IsLower ? gene1 : gene2;
            GeneRangeData upperGene = region1IsLower ? gene2 : gene1;

            GenePhaseRegion lowerRegion = region1IsLower ? region1 :region2;
            GenePhaseRegion upperRegion = region1IsLower ? region2 :region1;

            boolean lowerGeneIsUpstream = (isDel == isForwardStrand);
            boolean upperGeneIsUpstream = !lowerGeneIsUpstream;

            boolean phaseMatched = hasAnyPhaseMatch(lowerRegion, upperRegion, false);

            if (!phaseMatched)
            {
                if (upperGeneIsUpstream && regionsPhaseMatched(upperRegion, lowerRegion))
                    phaseMatched = true;
                else if (!upperGeneIsUpstream && regionsPhaseMatched(lowerRegion, upperRegion))
                    phaseMatched = true;
            }

            if (!phaseMatched)
                continue;

            foundMatch = true;

            Map<Integer, Long> bucketOverlapCounts = calcOverlapBucketAreas(
                    mProximateBucketLengths, lowerGene, upperGene, lowerRegion, upperRegion, isDel);

            if (!bucketOverlapCounts.isEmpty())
            {
                if (mLogVerbose)
                {
                    LOGGER.debug("gene({}: {}) and gene({}: {}) matched-region lower({} -> {}) upper({} -> {}) phase({}):",
                            lowerGene.GeneData.GeneId, lowerGene.GeneData.GeneName,
                            upperGene.GeneData.GeneId, upperGene.GeneData.GeneName,
                            lowerRegion.start(), lowerRegion.end(), upperRegion.start(), upperRegion.end(),
                            lowerRegion.getCombinedPhase());
                }

                for (Map.Entry<Integer, Long> entry : bucketOverlapCounts.entrySet())
                {
                    int bucketIndex = entry.getKey();
                    long overlap = entry.getValue();
                    addGeneFusionData(lowerGene, upperGene, overlap, isDel, bucketIndex);

                    if (mLogVerbose)
                    {
                        long bucketLen = mProximateBucketLengths.get(bucketIndex);
                        LOGGER.debug("matched-region bucketLength({}: {}) overlap({})", bucketIndex, bucketLen, overlap);
                    }
                }
            }
        }

        return foundMatch;
    }


    public static final int BUCKET_MIN = 0;
    public static final int BUCKET_MAX = 1;

    public long[] getBucketLengthMinMax(boolean isDel, int bucketIndex)
    {
        if(bucketIndex >= mProximateBucketLengths.size())
            return new long[2];

        return new long[] {mProximateBucketLengths.get(bucketIndex), mProximateBucketLengths.get(bucketIndex + 1)};
    }

    public void addGeneFusionData(final GeneRangeData lowerGene, final GeneRangeData upperGene, long overlapCount, boolean isDel, int bucketIndex)
    {
        long[] bucketMinMax = getBucketLengthMinMax(isDel, bucketIndex);
        long bucketWidth = bucketMinMax[BUCKET_MAX] - bucketMinMax[BUCKET_MIN];
        double fusionRate = overlapCount / (bucketWidth * GENOME_BASE_COUNT);

        if(fusionRate < MIN_FUSION_RATE)
            return;

        if(mLogVerbose)
        {
            LOGGER.debug("gene pair({} & {}) adding {} overlap({}) for bucket length index({})",
                    lowerGene.GeneData.GeneName, upperGene.GeneData.GeneName,
                    isDel ? "DEL" : "DUP", overlapCount, bucketIndex);
        }

        Map<String, Map<Integer,Long>> genePairCounts = isDel ? mDelGenePairCounts : mDupGenePairCounts;

        final String genePair = lowerGene.GeneData.GeneId + GENE_PAIR_DELIM + upperGene.GeneData.GeneId;
        Map<Integer,Long> bucketOverlapCounts = genePairCounts.get(genePair);

        if(bucketOverlapCounts == null)
        {
            bucketOverlapCounts = Maps.newHashMap();
            genePairCounts.put(genePair, bucketOverlapCounts);
        }

        setBucketLengthData(bucketOverlapCounts, bucketIndex, overlapCount);
    }

    public void logGlobalCounts()
    {
        LOGGER.info("GLOBAL_COUNTS: Type,GlobalCount,LengthMin,LengthMax");
        LOGGER.info("GLOBAL_COUNTS: LONG_DDI,{},{},{}", mGlobalLongDelDupInvCount, 5000000, 5000000);
        LOGGER.info("GLOBAL_COUNTS: INV,{},{},{}", mGlobalShortInvCount, 1000, 2000);

        // assign to correct bucket length
        for(int b = 0; b < mProximateBucketLengths.size() - 1; ++b)
        {
            long minLength = mProximateBucketLengths.get(b);
            long maxLength = mProximateBucketLengths.get(b + 1);

            LOGGER.info("GLOBAL_COUNTS: DEL,{},{},{}", mGlobalProximateCounts.get(b), minLength, maxLength);
            LOGGER.info("GLOBAL_COUNTS: DUP,{},{},{}", mGlobalProximateCounts.get(b), minLength, maxLength);
        }
    }


}
