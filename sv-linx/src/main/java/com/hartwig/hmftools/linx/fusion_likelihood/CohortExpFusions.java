package com.hartwig.hmftools.linx.fusion_likelihood;

import static java.lang.Math.abs;
import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.round;

import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_P;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_ARM_Q;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.CHROMOSOME_LENGTHS;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getChromosomalArmLength;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.makeChrArmStr;
import static com.hartwig.hmftools.linx.fusion.FusionFinder.TRANSCRIPT_PROTEIN_CODING;
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
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.overlapsOtherRegions;
import static com.hartwig.hmftools.linx.fusion_likelihood.PhaseRegionUtils.splitOverlappingPhaseRegion;
import static com.hartwig.hmftools.linx.fusion_likelihood.RegionAllocator.DEFAULT_REGION_GRID_SIZE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.structural.annotation.EnsemblGeneData;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;
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

    private boolean mLogVerbose;

    public static int PRE_GENE_3P_DISTANCE = 10000;
    public static int SHORT_INV_BUCKET = 100000;
    public static int MIN_BUCKET_LENGTH = 100;
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
        mGlobalProximateCounts = Lists.newArrayList();
        mGlobalShortInvCount = 0;
        mGlobalLongDelDupInvCount = 0;
        mArmLengthFactor = 0;
        mLogVerbose = false;
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
                final List<TranscriptData> transDataList = geneTransCache.getTranscripts(geneData.GeneId);

                if (transDataList == null)
                    continue;

                generatePhaseRegions(geneRangeData, transDataList, geneTransCache);

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
                }
            }

            if(!armGeneList.isEmpty())
            {
                processArmGenes(chromosome, currentArm, armGeneList, armGeneEndFirstList);
            }

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
        divideOverlappingRegions(chromosome, arm, armGeneList);

        for(int i = 0; i <= 0; ++i)
        {
            int strand = (i == 0) ? 1 : -1;

            List<GeneRangeData> geneDataList = strand == 1 ? armGeneList : armGeneEndFirstList;

            LOGGER.info("chr({}) arm({}) finding same-gene fusions", chromosome, arm);

            for(GeneRangeData geneData : geneDataList)
            {
                if(geneData.GeneData.Strand != strand)
                    continue;

                if(geneData.hasProteinCoding())
                {
                    generateSameGeneCounts(geneData);
                }
            }

            LOGGER.info("chr({}) arm({}) finding proximate fusions", chromosome, arm);

            generateProximateCounts(geneDataList, strand);

            LOGGER.info("chr({}) arm({}) finding non-proximate fusions", chromosome, arm);

            generateNonProximateCounts(geneDataList, strand);
        }

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
                final List<TranscriptData> transDataList = geneTransCache.getTranscripts(geneData.GeneId);

                if (transDataList == null || transDataList.isEmpty())
                    continue;

                /*
                if(geneData.GeneId.equals("ENSG00000260411"))
                {
                    LOGGER.info("spec gene({})", geneData.GeneId);
                }
                */

                generatePhaseRegions(geneRangeData, transDataList, geneTransCache);

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
                }
            }
        }
    }

    public void generatePhaseRegions(
            GeneRangeData geneRangeData, final List<TranscriptData> transDataList, final SvGeneTranscriptCollection geneTransCache)
    {
        // convert each transcript's exons into a set of phase regions spanning each intronic section
        // also take each transcript and look for potential same-gene fusions
        List<GenePhaseRegion> phaseRegions = Lists.newArrayList();
        List<GenePhaseRegion> intronicPhaseRegions = Lists.newArrayList();
        List<GenePhaseRegion> allTranscriptRegions = Lists.newArrayList();

        for(TranscriptData transcript : transDataList)
        {
            int transId = transcript.TransId;
            long precedingGeneSAPos = geneTransCache.findPrecedingGeneSpliceAcceptorPosition(transId);

            List<GenePhaseRegion> transcriptRegions = createPhaseRegionsFromTranscript(geneRangeData.GeneData, transcript, precedingGeneSAPos);

            transcriptRegions.forEach(x -> x.setTransId(transId));

            // add the new transcript's set of regions only if they aren't a very close overlap with any existing regions of the same phase
            for(GenePhaseRegion transRegion : transcriptRegions)
            {
                if(!overlapsOtherRegions(transRegion, allTranscriptRegions, true, 0.75))
                    allTranscriptRegions.add(transRegion);
            }

            // only looking at canonical simplies the overlap logic but misses about 1/2 the actual same-gene fusions seen in prod
            // if(transcript.IsCanonical)
            //    geneRangeData.setTranscriptPhaseRegions(transcriptRegions);

            // consolidate regions where phases and pre-gene status overlap
            transcriptRegions.stream().forEach(x -> checkAddCombinedGenePhaseRegion(x, phaseRegions));
            // transcriptRegions.stream().forEach(x -> combineGeneIntronicPhaseRegion(x, intronicPhaseRegions));
        }

        geneRangeData.setTranscriptPhaseRegions(allTranscriptRegions);

        mergePhaseRegions(phaseRegions);
        geneRangeData.setPhaseRegions(phaseRegions);
        // geneRangeData.setIntronicPhaseRegions(intronicPhaseRegions);
    }

    public static List<GenePhaseRegion> createPhaseRegionsFromTranscript(
            final EnsemblGeneData geneData, final TranscriptData transcript, long precSpliceAcceptorPos)
    {
        // 5-prime rules - must be post-promotor
        // 3-prime rules: must be coding and > 1 exon, needs to find first splice acceptor and then uses its phasing
        // need to mark coding and non-coding regions

        List<GenePhaseRegion> transcriptRegions = Lists.newArrayList();

        // skip single-exon transcripts since without an intronic section their fusion likelihood is negligible
        if(transcript.exons().size() == 1)
            return transcriptRegions;

        boolean proteinCoding = transcript.BioType.equals(TRANSCRIPT_PROTEIN_CODING);

        for (int i = 0; i < transcript.exons().size() - 1; ++i)
        {
            ExonData exonData = transcript.exons().get(i);

            if (transcript.CodingStart == null)
            {
                // mark the whole transcript as a single UTR
                GenePhaseRegion phaseRegion = new GenePhaseRegion(
                        geneData.GeneId, transcript.TransStart, transcript.TransEnd, PHASE_NON_CODING);

                transcriptRegions.add(phaseRegion);
                break;
            }

            ExonData nextExonData = transcript.exons().get(i+1);

            if(geneData.Strand == 1 && nextExonData.ExonStart > transcript.CodingEnd)
                break;
            else if(geneData.Strand == -1 && exonData.ExonEnd < transcript.CodingStart)
                continue;

            // add an upstream gene region (only applicable for downstream fusion genes)
            // with a phase of -1 unless coding starts in the first exon or on the first base of the second exon
            if(precSpliceAcceptorPos > 0 && geneData.Strand == 1 && i == 0)
            {
                int regionPhase = (transcript.CodingStart < exonData.ExonEnd) ? nextExonData.ExonPhase : -1;

                long preDistance = max(transcript.TransStart - precSpliceAcceptorPos, 0);
                long upstreamDistance = min(preDistance, PRE_GENE_3P_DISTANCE);
                long regionStart = transcript.TransStart - upstreamDistance;
                long regionEnd = exonData.ExonStart - 1;

                if(regionStart < regionEnd)
                {
                    GenePhaseRegion phaseRegion =
                            new GenePhaseRegion(geneData.GeneId, regionStart, regionEnd, mapExonPhase(regionPhase));
                    phaseRegion.setPreGene(true, phaseRegion.Phase);

                    transcriptRegions.add(phaseRegion);
                }
            }
            else if(precSpliceAcceptorPos > 0 && geneData.Strand == -1 && i == transcript.exons().size() - 2)
            {
                int regionPhase = (transcript.CodingEnd > nextExonData.ExonStart) ? exonData.ExonPhase : -1;

                long regionStart = nextExonData.ExonEnd + 1;

                long preDistance = max(precSpliceAcceptorPos - transcript.TransEnd, 0);
                long upstreamDistance = min(preDistance, PRE_GENE_3P_DISTANCE);
                long regionEnd = transcript.TransEnd + upstreamDistance;

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

        transcriptRegions.forEach(x -> x.setProteinCoding(proteinCoding));

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

        regions.add(GenePhaseRegion.from(newRegion, 0, newRegion.length()));
    }

    public void generateSameGeneCounts(GeneRangeData geneData)
    {
        // look for skipped or repeated exons within the same transcript and within phased regions
        // use the region allocator to avoid double-counting the same 2 regions across different transcripts
        final List<GenePhaseRegion> transcriptRegions = geneData.getTranscriptPhaseRegions();

        if(transcriptRegions.isEmpty())
            return;

        // allow for
        // long geneLength = geneData.GeneData.GeneEnd - geneData.GeneData.GeneStart;
        // double rawBlockSize = pow(10, round(log10(geneLength))) * 0.01;
        // int blockSize = DEFAULT_REGION_GRID_SIZE; // max(DEFAULT_REGION_GRID_SIZE, (int)rawBlockSize);

        int bucketLengths = mProximateBucketLengths.size() - 1;
        RegionAllocator[] regionAllocators = new RegionAllocator[bucketLengths];

        for(int i = 0; i < bucketLengths; ++i)
        {
            int blockSize = (int)(mProximateBucketLengths.get(i) / 10);
            blockSize = max(blockSize, MIN_BUCKET_LENGTH);
            regionAllocators[i] = new RegionAllocator(blockSize);
        }

        for (int i = 0; i < transcriptRegions.size(); ++i)
        {
            GenePhaseRegion region1 = transcriptRegions.get(i);

            if (!region1.hasPhasedType())
                continue;

            for (int j = i + 1; j < transcriptRegions.size(); ++j)
            {
                GenePhaseRegion region2 = transcriptRegions.get(j);

                if (region1.transId() != region2.transId())
                    break;

                if (!region2.hasPhasedType())
                    continue;

                testProximatePhaseRegions(geneData, geneData, region1, region2, regionAllocators);
            }
        }
    }

    private void divideOverlappingRegions(final String chromosome, final String arm, List<GeneRangeData> geneRangeList)
    {
        int phaseRegions = geneRangeList.stream().mapToInt(x -> x.getPhaseRegions().size()).sum();
        int regionsRemoved = 0;

        // find all regions with an overlap, to later control their phase-match allocation
        for (int lgIndex = 0; lgIndex < geneRangeList.size(); ++lgIndex)
        {
            GeneRangeData lowerGene = geneRangeList.get(lgIndex);

            List<GenePhaseRegion> lowerRegions = lowerGene.getPhaseRegions();

            // don't allow same-gene fusions (they are handled within a transcript), so start the index at the next gene
            for (int ugIndex = lgIndex + 1; ugIndex < geneRangeList.size(); ++ugIndex)
            {
                GeneRangeData upperGene = geneRangeList.get(ugIndex);

                if (upperGene.GeneData.Strand != lowerGene.GeneData.Strand)
                    continue;

                if (!upperGene.Arm.equals(lowerGene.Arm))
                    break;

                /*
                if(lowerGene.GeneData.GeneId.equals("ENSG00000121236") && upperGene.GeneData.GeneId.equals("ENSG00000258588"))
                {
                    LOGGER.info("spec genes");
                }
                */

                List<GenePhaseRegion> upperRegions = upperGene.getPhaseRegions();

                int lrIndex = 0;

                while(lrIndex < lowerRegions.size())
                {
                    GenePhaseRegion lowerRegion = lowerRegions.get(lrIndex);

                    int lrCount = lowerRegions.size();

                    int urIndex = 0;

                    while(urIndex < upperRegions.size())
                    {
                        GenePhaseRegion upperRegion = upperRegions.get(urIndex);

                        if (!haveOverlap(lowerRegion, upperRegion, -PERMITTED_REGION_OVERLAP))
                        {
                            ++urIndex;
                            continue;
                        }

                        int urCount = upperRegions.size();

                        if(lrIndex >= lowerRegions.size() || urIndex >= upperRegions.size())
                        {
                            LOGGER.error("genes({} & {}) index errors", lowerGene.GeneData.GeneId, upperGene.GeneData.GeneId);
                            return;
                        }

                        splitOverlappingPhaseRegion(lowerRegion, lrIndex, lowerRegions, upperRegion, urIndex, upperRegions);

                        // check for a region removed from either list
                        if(urCount == upperRegions.size())
                            ++urIndex;
                        else
                            ++regionsRemoved;

                        if(lowerRegions.size() < lrCount)
                        {
                            ++regionsRemoved;
                            break;
                        }
                    }

                    if(lrCount == lowerRegions.size())
                        ++lrIndex;
                }
            }
        }

        int newPhaseRegions = geneRangeList.stream().mapToInt(x -> x.getPhaseRegions().size()).sum();

        int added = max(newPhaseRegions - phaseRegions + regionsRemoved, 0);

        LOGGER.debug("chromosome({}) arm({}) dividing phase regions: initial({}) removed({}) added({}) final({})",
                chromosome, arm, phaseRegions, regionsRemoved, added, newPhaseRegions);
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
                    else if(maxDistance <= proximateLimit)
                    {
                        type = NON_PROX_TYPE_MEDIUM_INV;
                    }
                    else
                    {
                        type = NON_PROX_TYPE_LONG_SAME_ARM;
                    }
                }

                // calculate phasing overlap areas
                for (GenePhaseRegion region1 : gene1.getPhaseRegions())
                {
                    // the downstream gene of the potential fusion cannot be non-coding
                    for (GenePhaseRegion region2 : gene2.getPhaseRegions())
                    {
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

                        long regionOverlap = region1.length() * region2.length();

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
        // first tally up all phase counts per chromosome to then be used against all the others

        // cache the total counts (expressed as a length) for each unique phase combination per chromosome
        Map<String, List<GenePhaseRegion>> chrPhaseRegionsPosStrand = Maps.newHashMap();
        Map<String, List<GenePhaseRegion>> chrPhaseRegionsNegStrand = Maps.newHashMap();

        for(Map.Entry<String, List<GeneRangeData>> entry : mChrGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            final List<GeneRangeData> geneList = entry.getValue();

            if(geneList.isEmpty())
                continue;

            List<GenePhaseRegion> phaseRegionsPosStrand = Lists.newArrayList();
            List<GenePhaseRegion> phaseRegionsNegStrand = Lists.newArrayList();

            for(final GeneRangeData geneData : geneList)
            {
                for(GenePhaseRegion region : geneData.getPhaseRegions())
                {
                    // add to the chromosome's summary phasing regions by strand, without concern for overlaps
                    // across genes on the same strand
                    if(geneData.GeneData.Strand == 1)
                    {
                        addPhaseRegion(phaseRegionsPosStrand, region);
                    }
                    else
                    {
                        addPhaseRegion(phaseRegionsNegStrand, region);
                    }
                }
            }

            chrPhaseRegionsPosStrand.put(chromosome, phaseRegionsPosStrand);
            chrPhaseRegionsNegStrand.put(chromosome, phaseRegionsNegStrand);
        }

        // now test each chromosome's genes against the sum across all others
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
                            (j == 0) ? chrPhaseRegionsPosStrand : chrPhaseRegionsNegStrand;

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

    public static int PERMITTED_REGION_OVERLAP = 5;

    public boolean testProximatePhaseRegions(GeneRangeData gene1, GeneRangeData gene2, GenePhaseRegion region1, GenePhaseRegion region2)
    {
        return testProximatePhaseRegions(gene1, gene2, region1, region2, null);
    }

    public boolean testProximatePhaseRegions(GeneRangeData gene1, GeneRangeData gene2, GenePhaseRegion region1, GenePhaseRegion region2,
            RegionAllocator[] regionAllocators)
    {
        // ignore overlapping regions for now since it's not clear whether a DUP or DEL would be required
        if (haveOverlap(region1, region2, -PERMITTED_REGION_OVERLAP)) // allow bases with an exact base overlap through
            return false;

        // skip tiny overlaps
        if(region1.length() < PERMITTED_REGION_OVERLAP || region2.length() < PERMITTED_REGION_OVERLAP)
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

            Map<Integer, Long> bucketOverlapCounts = calcOverlapBucketAreas(
                    mProximateBucketLengths, regionAllocators, lowerGene, upperGene, lowerRegion, upperRegion, isDel);

            if (!bucketOverlapCounts.isEmpty())
            {
                foundMatch = true;

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
