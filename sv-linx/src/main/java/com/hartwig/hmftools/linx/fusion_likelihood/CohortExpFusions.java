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
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.hasAnyPhaseMatch;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.haveOverlap;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.mapExonPhase;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseRegion.regionsPhaseMatched;
import static com.hartwig.hmftools.linx.fusion_likelihood.GenePhaseType.PHASE_NON_CODING;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_LONG_SAME_ARM;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_MEDIUM_INV;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_REMOTE;
import static com.hartwig.hmftools.linx.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_SHORT_INV;
import static com.hartwig.hmftools.linx.fusion_likelihood.LikelihoodCalc.calcOverlapBucketAreas;
import static com.hartwig.hmftools.linx.fusion_likelihood.LikelihoodCalc.setBucketLengthData;
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.checkAddCombinedGenePhaseRegion;
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.mergePhaseRegions;
import static com.hartwig.hmftools.linx.fusion_likelihood.RegionAllocator.DEFAULT_REGION_GRID_SIZE;
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

    private final Map<String, List<GeneRangeData>> mChrGeneDataMap;
    private Map<String, GeneRangeData> mGeneIdRangeDataMap;

    private Map<String,Map<Integer,Long>> mDelGenePairCounts; // pair of gene-pairs to their overlap counts keyed by bucket index
    private Map<String,Map<Integer,Long>> mDupGenePairCounts;

    // global counts by type and buck length
    private List<Integer> mGlobalProximateCounts; // indexed as per the proximate lengths
    private int mGlobalShortInvCount;
    private int mGlobalLongDelDupInvCount;
    private long mArmLengthFactor;

    private RegionAllocator mRegionAllocator;

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
        mChrGeneDataMap = Maps.newHashMap();
        mGeneIdRangeDataMap = null;
        mDelGenePairCounts = Maps.newHashMap();
        mDupGenePairCounts = Maps.newHashMap();
        mChrPhaseRegionsPosStrand = Maps.newHashMap();
        mChrPhaseRegionsNegStrand = Maps.newHashMap();
        mGlobalProximateCounts = Lists.newArrayList();
        mGlobalShortInvCount = 0;
        mGlobalLongDelDupInvCount = 0;
        mArmLengthFactor = 0;
        mLogVerbose = false;

        mRegionAllocator = new RegionAllocator(DEFAULT_REGION_GRID_SIZE);
    }

    public final Map<String, List<GeneRangeData>> getChrGeneRangeDataMap() { return mChrGeneDataMap; }
    public Map<String,Map<Integer,Long>> getDelGenePairCounts() { return mDelGenePairCounts; }
    public Map<String,Map<Integer,Long>> getDupGenePairCounts() { return mDupGenePairCounts; }

    public long getArmLengthFactor() { return mArmLengthFactor; }

    public void initialise(List<Long> proximateBucketLengths, int shortInvBucketLength)
    {
        mProximateBucketLengths = proximateBucketLengths;

        mProximateBucketLengths.stream().forEach(x -> mGlobalProximateCounts.add(0));

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
                mArmLengthFactor += armLength * (armLength - maxBucketLength);
        }

        LOGGER.info("arm length factor: {}", mArmLengthFactor);
    }

    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }

    public void initialiseGeneIdRangeDataMap()
    {
        mGeneIdRangeDataMap = Maps.newHashMap();
    }

    public GeneRangeData findGeneRangeData(final String geneId)
    {
        return mGeneIdRangeDataMap != null ? mGeneIdRangeDataMap.get(geneId) : null;
    }

    public void generateExpectedFusions(final SvGeneTranscriptCollection geneTransCache,
            List<String> restrictedChromosomes, List<String> restrictedGeneIds)
    {
        /* for each chromosome and arm:
            - convert each gene's transcripts into a set of phase regions
            - mark any overlapping phase regions between genes
            - determine same-gene fusions
            - determine proximate fusions
            - determine non-proximate fusions
            - finally determine global fusion
        */

        final Map<String, List<EnsemblGeneData>> chrGeneDataMap = geneTransCache.getChrGeneDataMap();

        for(Map.Entry<String, List<EnsemblGeneData>> entry : chrGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            if(!restrictedChromosomes.isEmpty() && !restrictedChromosomes.contains(chromosome))
                continue;

            final List<EnsemblGeneData> chrGenes = entry.getValue();

            if(chrGenes.isEmpty())
                continue;

            LOGGER.info("generating phase counts for chromosome({})", chromosome);

            List<GeneRangeData> chrGeneList = Lists.newArrayList();
            List<GeneRangeData> armGeneList = Lists.newArrayList();
            List<GeneRangeData> armGeneEndFirstList = Lists.newArrayList();
            List<GenePhaseRegion> chrPhaseRegionsPosStrand = Lists.newArrayList();
            List<GenePhaseRegion> chrPhaseRegionsNegStrand = Lists.newArrayList();

            String currentArm = "";

            for(final EnsemblGeneData geneData :chrGenes)
            {
                if(!restrictedGeneIds.isEmpty() && !restrictedGeneIds.contains(geneData.GeneId))
                    continue;

                GeneRangeData geneRangeData = new GeneRangeData(geneData);

                if(currentArm == "")
                {
                    currentArm = geneRangeData.Arm;
                }
                else if(!currentArm.equals(geneRangeData.Arm))
                {
                    processArmGenes(chromosome, currentArm, armGeneList, armGeneEndFirstList);
                    armGeneList.clear();
                    armGeneEndFirstList.clear();

                    currentArm = geneRangeData.Arm;
                }

                chrGeneList.add(geneRangeData);
                armGeneList.add(geneRangeData);

                int index = 0;
                for(; index < armGeneEndFirstList.size(); ++index)
                {
                    final GeneRangeData rgd = armGeneEndFirstList.get(index);

                    if(geneData.GeneEnd < rgd.GeneData.GeneEnd)
                        break;
                }

                armGeneEndFirstList.add(index, geneRangeData);

                // load from Ensembl transcript and exon data
                final List<TranscriptExonData> transExonDataList = geneTransCache.getTransExonData(geneData.GeneId);

                if (transExonDataList == null)
                    continue;

                generatePhaseRegions(geneRangeData, transExonDataList, geneTransCache);

                for(GenePhaseRegion region : geneRangeData.getPhaseRegions())
                {
                    if(mLogVerbose)
                    {
                        LOGGER.debug("gene({}) {}", region.GeneId, region.toString());
                    }

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
                        addPhaseRegion(chrPhaseRegionsPosStrand, region);
                    }
                    else
                    {
                        addPhaseRegion(chrPhaseRegionsNegStrand, region);
                    }
                }
            }

            if(!armGeneList.isEmpty())
            {
                processArmGenes(chromosome, currentArm, armGeneList, armGeneEndFirstList);
            }

            mChrPhaseRegionsPosStrand.put(chromosome, chrPhaseRegionsPosStrand);
            mChrPhaseRegionsNegStrand.put(chromosome, chrPhaseRegionsNegStrand);

            mChrGeneDataMap.put(chromosome, chrGeneList);
        }

        // finally generate remote / translocation counts
        generateRemoteCounts();
    }

    private void processArmGenes(final String chromosome, final String arm,
            List<GeneRangeData> armGeneList, List<GeneRangeData> armGeneEndFirstList)
    {
        LOGGER.info("chr({}) arm({}) processing {} genes", chromosome, arm, armGeneList.size());

        // mark all overlapping regions to avoid reallocation of base overlaps
        markOverlappingRegions(chromosome, arm, armGeneList);

        for(int i = 0; i <= 0; ++i)
        {
            int strand = (i == 0) ? 1 : -1;

            // reset the allocation tracker since the scope of its state is per arm and strand
            if(mRegionAllocator.allocationCount() > 0)
            {
                LOGGER.info("chr({}) arm({}) clearing {} region allocations",
                        chromosome, arm, armGeneList.size(), mRegionAllocator.allocationCount());
            }

            mRegionAllocator.reset();

            LOGGER.info("chr({}) arm({}) finding same-gene fusions", chromosome, arm);

            for(GeneRangeData geneData : armGeneList)
            {
                if(geneData.GeneData.Strand != strand)
                    continue;

                generateSameGeneCounts(geneData);
            }

            LOGGER.info("chr({}) arm({}) finding proximate fusions", chromosome, arm);

            generateProximateCounts(armGeneList, strand);

            LOGGER.info("chr({}) arm({}) finding non-proximate fusions", chromosome, arm);

            generateNonProximateCounts(armGeneList, strand);
        }

        // check for INVs across both strands
        if(mRegionAllocator.allocationCount() > 0)
        {
            LOGGER.info("chr({}) arm({}) clearing {} region allocations",
                    chromosome, arm, armGeneList.size(), mRegionAllocator.allocationCount());
        }

        mRegionAllocator.reset();

        generateNonProximateCounts(armGeneList, 0);
    }

    public void generateGenePhasingCounts(final SvGeneTranscriptCollection geneTransCache,
            List<String> restrictedChromosomes, List<String> restrictedGeneIds)
    {
        // for each gene, walk through all its transcripts and count up the number of bases in each phasing region (eg 5'UTR, 0-2)
        // also convert these into combined phasing regions (eg where say phase 0 and 1 overlap) for using in downstream analysis
        // finally gather up the total bases by phasing type across each chromosome
        final Map<String, List<EnsemblGeneData>> chrGeneDataMap = geneTransCache.getChrGeneDataMap();

        for(Map.Entry<String, List<EnsemblGeneData>> entry : chrGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            if(!restrictedChromosomes.isEmpty() && !restrictedChromosomes.contains(chromosome))
                continue;

            LOGGER.info("generating phase counts for chromosome({})", chromosome);

            List<GeneRangeData> geneList = Lists.newArrayList();
            List<GeneRangeData> geneEndFirstList = Lists.newArrayList();
            List<GenePhaseRegion> chrPhaseRegionsPosStrand = Lists.newArrayList();
            List<GenePhaseRegion> chrPhaseRegionsNegStrand = Lists.newArrayList();
            long posStrandTotal = 0;
            long negStrandTotal = 0;

            for(final EnsemblGeneData geneData :entry.getValue())
            {
                if(!restrictedGeneIds.isEmpty() && !restrictedGeneIds.contains(geneData.GeneId))
                    continue;

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

                // load from Ensembl transcript and exon data
                final List<TranscriptExonData> transExonDataList = geneTransCache.getTransExonData(geneData.GeneId);

                if (transExonDataList == null)
                    continue;

                /*
                if(geneData.GeneId.equals("ENSG00000260411"))
                {
                    LOGGER.info("spec gene({})", geneData.GeneId);
                }
                */

                generatePhaseRegions(geneRangeData, transExonDataList, geneTransCache);

                for(GenePhaseRegion region : geneRangeData.getPhaseRegions())
                {
                    if(mLogVerbose)
                    {
                        LOGGER.debug("gene({}) {}", region.GeneId, region.toString());
                    }

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

            // markOverlappingRegions(geneList);

            mChrPhaseRegionsPosStrand.put(chromosome, chrPhaseRegionsPosStrand);
            mChrPhaseRegionsNegStrand.put(chromosome, chrPhaseRegionsNegStrand);

            LOGGER.debug("chr({}) posStrandTotal({}) negStrandTotal({})", chromosome, posStrandTotal, negStrandTotal);

            mChrGeneDataMap.put(chromosome, geneList);
            // mChrReverseGeneDataMap.put(chromosome, geneEndFirstList);
        }
    }

    public void generatePhaseRegions(
            GeneRangeData geneRangeData, final List<TranscriptExonData> transExonDataList, final SvGeneTranscriptCollection geneTransCache)
    {
        // convert each transcript's exons into a set of phase regions spanning each intronic section
        // also take each transcript and look for potential same-gene fusions
        int teIndex = 0;
        List<TranscriptExonData> transcriptExons = nextTranscriptExons(transExonDataList, teIndex);

        List<GenePhaseRegion> matchedTransRegions = Lists.newArrayList();
        List<GenePhaseRegion> phaseRegions = Lists.newArrayList();
        List<GenePhaseRegion> intronicPhaseRegions = Lists.newArrayList();

        while (!transcriptExons.isEmpty())
        {
            int transId = transcriptExons.get(0).TransId;
            long precedingGeneSAPos = geneTransCache.findPrecedingGeneSpliceAcceptorPosition(transId);

            List<GenePhaseRegion> transcriptRegions = createPhaseRegionsFromTranscript(geneRangeData.GeneData, transcriptExons, precedingGeneSAPos);

            transcriptRegions.forEach(x -> x.setTransId(transId));

            geneRangeData.setTranscriptPhaseRegions(transcriptRegions);

            // consolidate regions where phases and pre-gene status overlap
            transcriptRegions.stream().forEach(x -> checkAddCombinedGenePhaseRegion(x, phaseRegions));
            // transcriptRegions.stream().forEach(x -> combineGeneIntronicPhaseRegion(x, intronicPhaseRegions));

            // generateSameGeneCounts(geneRangeData, transcriptRegions);
            // generateSameGeneCounts(geneRangeData, transcriptRegions, matchedTransRegions);

            teIndex += transcriptExons.size();
            transcriptExons = nextTranscriptExons(transExonDataList, teIndex);
        }

        mergePhaseRegions(phaseRegions);
        geneRangeData.setPhaseRegions(phaseRegions);
        // geneRangeData.setIntronicPhaseRegions(intronicPhaseRegions);
    }

    public static List<GenePhaseRegion> createPhaseRegionsFromTranscript(
            final EnsemblGeneData geneData, final List<TranscriptExonData> transcriptExons, long precSpliceAcceptorPos)
    {
        // 5-prime rules - must be post-promotor
        // 3-prime rules: must be coding and > 1 exon, needs to find first splice acceptor and then uses its phasing
        // need to mark coding and non-coding regions

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

                transcriptRegions.add(phaseRegion);
                break;
            }

            TranscriptExonData nextExonData = transcriptExons.get(i+1);

            if(geneData.Strand == 1 && nextExonData.ExonStart > exonData.CodingEnd)
                break;
            else if(geneData.Strand == -1 && exonData.ExonEnd < exonData.CodingStart)
                continue;

            // add an upstream gene region (only applicable for downstream fusion genes)
            // with a phase of -1 unless coding starts in the first exon or on the first base of the second exon
            if(precSpliceAcceptorPos > 0 && geneData.Strand == 1 && i == 0)
            {
                int regionPhase = (exonData.CodingStart < exonData.ExonEnd) ? nextExonData.ExonPhase : -1;

                long preDistance = max(exonData.TransStart - precSpliceAcceptorPos, 0);
                long upstreamDistance = min(preDistance, PRE_GENE_3P_DISTANCE);
                long regionStart = exonData.TransStart - upstreamDistance;
                long regionEnd = exonData.ExonStart - 1;

                if(regionStart < regionEnd)
                {
                    GenePhaseRegion phaseRegion =
                            new GenePhaseRegion(geneData.GeneId, regionStart, regionEnd, mapExonPhase(regionPhase));
                    phaseRegion.setPreGene(true, phaseRegion.Phase);

                    transcriptRegions.add(phaseRegion);
                }
            }
            else if(precSpliceAcceptorPos > 0 && geneData.Strand == -1 && i == transcriptExons.size() - 2)
            {
                int regionPhase = (exonData.CodingEnd > nextExonData.ExonStart) ? exonData.ExonPhase : -1;

                long regionStart = nextExonData.ExonEnd + 1;

                long preDistance = max(precSpliceAcceptorPos - exonData.TransEnd, 0);
                long upstreamDistance = min(preDistance, PRE_GENE_3P_DISTANCE);
                long regionEnd = nextExonData.TransEnd + upstreamDistance;

                if(regionStart < regionEnd)
                {
                    GenePhaseRegion phaseRegion =
                            new GenePhaseRegion(geneData.GeneId, regionStart, regionEnd, mapExonPhase(regionPhase));
                    phaseRegion.setPreGene(true, phaseRegion.Phase);

                    transcriptRegions.add(phaseRegion);
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

            transcriptRegions.add(phaseRegion);
        }

        return transcriptRegions;
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

    public void generateSameGeneCounts(GeneRangeData geneData)
    {
        final List<GenePhaseRegion> transcriptRegions = geneData.getTranscriptPhaseRegions();

        for (int i = 0; i < transcriptRegions.size(); ++i)
        {
            GenePhaseRegion region1 = transcriptRegions.get(i);

            for (int j = i+1; j < transcriptRegions.size(); ++j)
            {
                GenePhaseRegion region2 = transcriptRegions.get(j);

                if(region1.transId() != region2.transId())
                    break;

                testProximatePhaseRegions(geneData, geneData, region1, region2);
            }
        }
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

                for (final GenePhaseRegion region : matchedTransRegions)
                {
                    if (!haveOverlap(region, newRegion, 0))
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

    private void markOverlappingRegions(final String chromosome, final String arm, List<GeneRangeData> geneRangeList)
    {
        int overlaps = 0;

        // find all regions with an overlap, to later control their phase-match allocation
        for (int lowerIndex = 0; lowerIndex < geneRangeList.size(); ++lowerIndex)
        {
            GeneRangeData lowerGene = geneRangeList.get(lowerIndex);

            // don't allow same-gene fusions (they are handled within a transcript), so start the index at the next gene
            for (int upperIndex = lowerIndex + 1; upperIndex < geneRangeList.size(); ++upperIndex)
            {
                GeneRangeData upperGene = geneRangeList.get(upperIndex);

                if (upperGene.GeneData.Strand != lowerGene.GeneData.Strand)
                    continue;

                if (!upperGene.Arm.equals(lowerGene.Arm))
                    break;

                for (GenePhaseRegion lowerRegion : lowerGene.getPhaseRegions())
                {
                    if(lowerRegion.hasOverlaps())
                        continue;

                    for (GenePhaseRegion upperRegion : upperGene.getPhaseRegions())
                    {
                        if(upperRegion.hasOverlaps())
                            continue;

                        if (!haveOverlap(lowerRegion, upperRegion, 0))
                            continue;

                        lowerRegion.setHasOverlaps(true);
                        upperRegion.setHasOverlaps(true);
                        ++overlaps;
                    }
                }
            }
        }

        if(overlaps > 0)
        {
            LOGGER.info("chr({}) arm({}) has {} overlapping regions", chromosome, arm, overlaps);
        }
    }

    public long getMaxBucketLength()
    {
        return !mProximateBucketLengths.isEmpty() ? mProximateBucketLengths.get(mProximateBucketLengths.size() - 1) : MIN_BUCKET_LENGTH;
    }

    private void generateProximateCounts(List<GeneRangeData> geneList, int strandMatch)
    {
        long rangeLimit = getMaxBucketLength();

        for(int lowerIndex = 0; lowerIndex < geneList.size(); ++lowerIndex)
        {
            GeneRangeData lowerGene = geneList.get(lowerIndex);

            if(lowerGene.GeneData.Strand != strandMatch)
                continue;

            /*
            if(lowerGene.GeneData.GeneId.equals("ENSG00000258643"))
            {
                LOGGER.info("spec gene({})", lowerGene.GeneData.GeneId);
            }
            */

            // don't allow same-gene fusions (they are handled within a transcript), so start the index at the next gene
            for(int upperIndex = lowerIndex+1; upperIndex < geneList.size(); ++upperIndex)
            {
                GeneRangeData upperGene = geneList.get(upperIndex);

                if(upperGene.GeneData.Strand != strandMatch)
                    continue;

                // exit if the gene is now too far away from the lower gene
                if(upperGene.GeneData.GeneStart - lowerGene.GeneData.GeneEnd > rangeLimit)
                    break;

                if(!upperGene.Arm.equals(lowerGene.Arm))
                    break;

                // find all matching phasing regions and account for any duplicated overlapping regions from different phasing matches
                // for all those found, assign them to a length bucket and find the number of bases that fall into the length bucket region
                for (GenePhaseRegion lowerRegion : lowerGene.getPhaseRegions())
                {
                    for (GenePhaseRegion upperRegion : upperGene.getPhaseRegions())
                    {
                        testProximatePhaseRegions(lowerGene, upperGene, lowerRegion, upperRegion);

                        /*
                        if (!haveOverlap(lowerRegion, upperRegion, 0))
                        {
                            testProximatePhaseRegions(lowerGene, upperGene, lowerRegion, upperRegion);
                        }
                        else
                        {
                            testOverlappingProximateRegions(lowerGene, upperGene, lowerRegion, upperRegion);
                        }
                        */
                    }
                }
            }
        }
    }

    public void generateNonProximateCounts(List<GeneRangeData> geneList, int strandMatch)
    {
        // test for candidate fusions between long DELs and DUPs, and then for INVs of various lengths
        long proximateLimit = getMaxBucketLength();

        for(int i = 0; i < geneList.size(); ++i)
        {
            GeneRangeData gene1 = geneList.get(i);

            if(strandMatch != 0 && gene1.GeneData.Strand != strandMatch)
                continue;

            // and now the local ones outside the DEL and DUP proximity lengths, ignoring any strand checking
            // this is considering fusions between DUPs, DUPs or INVs
            for(int j = i + 1; j < geneList.size(); ++j)
            {
                GeneRangeData gene2 = geneList.get(j);

                if(!gene1.Arm.equals(gene2.Arm))
                    break;

                if(strandMatch == 0)
                {
                    // looking for INVs
                    if (gene1.GeneData.Strand == gene2.GeneData.Strand)
                        continue;
                }
                else
                {
                    // looking for DELs and DUPs
                    if (gene2.GeneData.Strand != strandMatch)
                        continue;
                }

                // find furthest distance between these 2 genes
                long maxDistance = max(abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneStart),
                        abs(gene1.GeneData.GeneEnd - gene2.GeneData.GeneEnd));

                maxDistance = max(maxDistance, abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneEnd));
                maxDistance = max(maxDistance, abs(gene1.GeneData.GeneEnd - gene2.GeneData.GeneStart));

                int type = -1;

                // assign the proximity to a category to build up global stats
                if (strandMatch != 0)
                {
                    if(maxDistance <= proximateLimit)
                    {
                        // proximate DELs and DUPs handled elsewhere
                        continue;
                    }
                    else
                    {
                        type = NON_PROX_TYPE_LONG_SAME_ARM;
                    }
                }
                else
                {
                    if(maxDistance < SHORT_INV_BUCKET)
                    {
                        type = NON_PROX_TYPE_SHORT_INV;
                    }
                    else
                    {
                        type = NON_PROX_TYPE_MEDIUM_INV;
                    }
                }

                // calculate phasing overlap areas
                for (GenePhaseRegion region1 : gene1.getPhaseRegions())
                {
                    // the downstream gene of the potential fusion cannot be non-coding
                    for (GenePhaseRegion region2 : gene2.getPhaseRegions())
                    {
                        boolean hasOverlaps = region1.hasOverlaps() && region2.hasOverlaps();

                        boolean matchGene1AsUp = false;
                        boolean matchGene1AsDown = false;

                        if(hasAnyPhaseMatch(region1, region2, false))
                        {
                            matchGene1AsUp = true;
                            matchGene1AsDown = true;
                        }
                        else
                        {
                            // check each direction in turn
                            if (regionsPhaseMatched(region1, region2))
                            {
                                matchGene1AsUp = true;
                            }

                            if (regionsPhaseMatched(region2, region1))
                            {
                                matchGene1AsDown = true;
                            }
                        }

                        if(!matchGene1AsUp && !matchGene1AsDown)
                            continue;

                        long regionOverlap = !hasOverlaps ? (region1.length() * region2.length()) :
                                mRegionAllocator.allocateBases(region1.start(), region1.end(), region2.start(), region2.end());

                        if(matchGene1AsUp)
                        {
                            gene1.addBaseOverlapCountUpstream(type, regionOverlap);
                            gene2.addBaseOverlapCountDownstream(type, regionOverlap);
                        }

                        if(matchGene1AsDown)
                        {
                            gene2.addBaseOverlapCountUpstream(type, regionOverlap);
                            gene1.addBaseOverlapCountDownstream(type, regionOverlap);
                        }
                    }
                }
            }
        }
    }

    private void generateRemoteCounts()
    {
        for (Map.Entry<String, List<GeneRangeData>> entry : mChrGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            LOGGER.info("calculating remote overlap counts for chr({})", chromosome);

            List<GeneRangeData> chrGeneList = entry.getValue();

            for(GeneRangeData gene : chrGeneList)
            {
                // determine the remote overlap count by comparing this gene's regions to all matching remote ones
                for (int j = 0; j <= 1; ++j)
                {
                    final Map<String, List<GenePhaseRegion>> chrPhaseRegions =
                            (j == 0) ? mChrPhaseRegionsPosStrand : mChrPhaseRegionsNegStrand;

                    for (Map.Entry<String, List<GenePhaseRegion>> chrEntry : chrPhaseRegions.entrySet())
                    {
                        if (chrEntry.getKey().equals(chromosome)) // skip the same chromosome
                            continue;

                        final List<GenePhaseRegion> remoteRegions = chrEntry.getValue();

                        for (GenePhaseRegion region : gene.getPhaseRegions())
                        {
                            for (GenePhaseRegion remoteRegion : remoteRegions)
                            {
                                // test gene as an downstream vs all remote upstream phasing regions
                                if (hasAnyPhaseMatch(remoteRegion, region, false) || regionsPhaseMatched(remoteRegion, region))
                                {
                                    gene.addBaseOverlapCountDownstream(NON_PROX_TYPE_REMOTE, region.length() * remoteRegion.length());
                                }

                                // then as an upstream partner
                                if (hasAnyPhaseMatch(region, remoteRegion, false) || regionsPhaseMatched(region, remoteRegion))
                                {
                                    gene.addBaseOverlapCountUpstream(NON_PROX_TYPE_REMOTE, region.length() * remoteRegion.length());
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
        if (haveOverlap(region1, region2, -1)) // allow bases with an exact base overlap through
            return false;

        // skip tiny overlaps
        if(region1.length() < 5 || region2.length() < 5)
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
                    mProximateBucketLengths, mRegionAllocator, lowerGene, upperGene, lowerRegion, upperRegion, isDel);

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

    public void testOverlappingProximateRegions(GeneRangeData lowerGene, GeneRangeData upperGene,
            GenePhaseRegion lowerRegion, GenePhaseRegion upperRegion)
    {
        // preliminary check for a phase match before going to effort of splitting regions
        if(!hasAnyPhaseMatch(lowerRegion, upperRegion, false) && !regionsPhaseMatched(lowerRegion, upperRegion)
                && !hasAnyPhaseMatch(upperRegion, lowerRegion,false) && !regionsPhaseMatched(upperRegion, lowerRegion))
        {
            return;
        }

        long overlapStart = max(lowerRegion.start(), upperRegion.start());
        long overlapEnd = min(lowerRegion.end(), upperRegion.end());

        // first process the regions before and after the overlapping region
        if(lowerRegion.start() < upperRegion.start())
        {
            GenePhaseRegion lowerSplit =
                    new GenePhaseRegion(lowerRegion.GeneId, lowerRegion.start(), upperRegion.start(),
                            lowerRegion.getPhaseArray(), lowerRegion.getPreGenePhaseStatus());

            testProximatePhaseRegions(lowerGene, upperGene, lowerSplit, upperRegion);
        }
        else if(upperRegion.start() < lowerRegion.start())
        {
            GenePhaseRegion lowerSplit =
                    new GenePhaseRegion(upperRegion.GeneId, upperRegion.start(), lowerRegion.start(),
                            upperRegion.getPhaseArray(), upperRegion.getPreGenePhaseStatus());

            testProximatePhaseRegions(upperGene, lowerGene, lowerSplit, lowerRegion);
        }

        if(lowerRegion.end() > upperRegion.end())
        {
            GenePhaseRegion upperSplit = new GenePhaseRegion(lowerRegion.GeneId, upperRegion.end(), lowerRegion.end(),
                    lowerRegion.getPhaseArray(), lowerRegion.getPreGenePhaseStatus());

            GenePhaseRegion lowerSplit =
                    new GenePhaseRegion(upperRegion.GeneId, overlapStart, overlapEnd,
                            upperRegion.getPhaseArray(), upperRegion.getPreGenePhaseStatus());

            testProximatePhaseRegions(upperGene, lowerGene, lowerSplit, upperSplit);
        }
        else if(upperRegion.end() > lowerRegion.end())
        {
            GenePhaseRegion upperSplit = new GenePhaseRegion(upperRegion.GeneId, lowerRegion.end(), upperRegion.end(),
                    upperRegion.getPhaseArray(), upperRegion.getPreGenePhaseStatus());

            GenePhaseRegion lowerSplit =
                    new GenePhaseRegion(lowerRegion.GeneId, overlapStart, overlapEnd,
                            lowerRegion.getPhaseArray(), lowerRegion.getPreGenePhaseStatus());

            testProximatePhaseRegions(lowerGene, upperGene, lowerSplit, upperSplit);
        }

        long minBucketLength = mProximateBucketLengths.get(0);

        // next split the overlap into smaller chunks to proces these one by one

        double segment = (overlapEnd - overlapStart) / 5.0;
        long segmentSize = max((long)floor(segment), minBucketLength);

        // break the regions up into non-overlapping pieces and process them one by one
        List<GenePhaseRegion> lowerSplits = Lists.newArrayList();
        List<GenePhaseRegion> upperSplits = Lists.newArrayList();

        for(int i = 0; i < 5; ++i)
        {
            long regionStart = overlapStart + (i * segmentSize);
            long regionEnd = regionStart + segmentSize;

            if(regionEnd > overlapEnd)
                break;

            lowerSplits.add(new GenePhaseRegion(lowerRegion.GeneId, regionStart, regionEnd,
                    lowerRegion.getPhaseArray(), lowerRegion.getPreGenePhaseStatus()));

            upperSplits.add(new GenePhaseRegion(upperRegion.GeneId, regionStart, regionEnd,
                    upperRegion.getPhaseArray(), upperRegion.getPreGenePhaseStatus()));
        }

        for(int i = 0; i < lowerSplits.size(); ++i)
        {
            GenePhaseRegion lowerSegment = lowerSplits.get(i);

            for(int j = 0; j < upperSplits.size(); ++j)
            {
                if(i == j)
                    continue;

                GenePhaseRegion upperSegment = upperSplits.get(j);
                testProximatePhaseRegions(lowerGene, upperGene, lowerSegment, upperSegment);
            }
        }
    }

    public void generateProximateFusionCounts()
    {
        /*
        if(mProximateBucketLengths.isEmpty())
            return;

        // walk through each chromosome, and test distances between each gene and the next ones for potential
        // fusions from DELs and DUPs within the configured bucket distances

        LOGGER.info("finding proximate fusion candidates on forward strand");

        for (Map.Entry<String, List<GeneRangeData>> entry : mChrGeneDataMap.entrySet())
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
        */
    }

    public void generateNonProximateCounts_old()
    {
        // for each arm and each gene in that arm sum up the overlapping counts against all genes beyond the specified
        // DEL and DUP max bucket length, and then all overlapping counts on remote arms
        long proximateLimit = getMaxBucketLength();

        for (Map.Entry<String, List<GeneRangeData>> entry : mChrGeneDataMap.entrySet())
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

                        for (GenePhaseRegion region : gene1.getPhaseRegions())
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

                    if(!gene1.ChromosomeArm.equals(gene2.ChromosomeArm))
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
                    if (j > i && maxDistance > proximateLimit)
                    {
                        ++mGlobalLongDelDupInvCount;
                        type = NON_PROX_TYPE_LONG_SAME_ARM;
                    }
                    else if(gene1.GeneData.Strand != gene2.GeneData.Strand)
                    {
                        if(maxDistance < SHORT_INV_BUCKET)
                        {
                            ++mGlobalShortInvCount;
                            type = NON_PROX_TYPE_SHORT_INV;
                        }
                        else
                        {
                            type = NON_PROX_TYPE_MEDIUM_INV;
                        }
                    }
                    else
                    {
                        // assign proximate genes to correct bucket length, same-genes will be counted here as well
                        for (int b = 0; b < mProximateBucketLengths.size() - 1; ++b)
                        {
                            long minLength = mProximateBucketLengths.get(b);
                            long maxLength = mProximateBucketLengths.get(b + 1);

                            if (minDistance > maxLength || maxDistance < minLength)
                                continue;

                            mGlobalProximateCounts.set(b, mGlobalProximateCounts.get(b) + 1);
                        }

                        continue;
                    }

                    // calculate phasing overlap areas
                    for (GenePhaseRegion region1 : gene1.getPhaseRegions())
                    {
                        // the downstream gene of the potential fusion cannot be non-coding
                        for (GenePhaseRegion region2 : gene2.getPhaseRegions())
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
