package com.hartwig.hmftools.svtools.fusion_likelihood;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;

import static com.hartwig.hmftools.common.gene.TranscriptProteinData.BIOTYPE_PROTEIN_CODING;
import static com.hartwig.hmftools.linx.analysis.SvUtilities.getChromosomalArmLength;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.P_ARM;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.Q_ARM;
import static com.hartwig.hmftools.common.purple.segment.ChromosomeArm.UNKNOWN;
import static com.hartwig.hmftools.svtools.fusion_likelihood.FusionLikelihood.FLC_LOGGER;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseRegion.hasAnyPhaseMatch;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseRegion.haveOverlap;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseRegion.mapExonPhase;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseRegion.regionsPhaseMatched;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GenePhaseType.PHASE_NON_CODING;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_LONG_SAME_ARM;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_MEDIUM_INV;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_REMOTE;
import static com.hartwig.hmftools.svtools.fusion_likelihood.GeneRangeData.NON_PROX_TYPE_SHORT_INV;
import static com.hartwig.hmftools.svtools.fusion_likelihood.LikelihoodCalc.calcOverlapBucketAreas;
import static com.hartwig.hmftools.svtools.fusion_likelihood.LikelihoodCalc.setBucketLengthData;
import static com.hartwig.hmftools.svtools.fusion_likelihood.PhaseRegionUtils.checkAddCombinedGenePhaseRegion;
import static com.hartwig.hmftools.svtools.fusion_likelihood.PhaseRegionUtils.divideOverlappingRegions;
import static com.hartwig.hmftools.svtools.fusion_likelihood.PhaseRegionUtils.mergePhaseRegions;
import static com.hartwig.hmftools.svtools.fusion_likelihood.PhaseRegionUtils.overlapsOtherRegions;
import static com.hartwig.hmftools.svtools.fusion_likelihood.RegionAllocator.DEFAULT_BUCKET_REGION_RATIO;
import static com.hartwig.hmftools.svtools.fusion_likelihood.RegionAllocator.MIN_BLOCK_SIZE;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.linx.analysis.SvUtilities;
import com.hartwig.hmftools.common.purple.segment.ChromosomeArm;

// routines for calculating expected fusion rates across the whole genome for a set of pre-described SV and length categories
public class CohortExpFusions
{
    // config
    private final List<Integer> mProximateBucketLengths;
    private int mMaxProximateBucketLength;

    private final Map<String, List<GeneRangeData>> mChrGeneDataMap;

    private final Map<String,Map<Integer,Long>> mDelGenePairCounts; // pair of gene-pairs to their overlap counts keyed by bucket index
    private final Map<String,Map<Integer,Long>> mDupGenePairCounts;

    // global counts by type and buck length
    private final List<Integer> mGlobalProximateCounts; // indexed as per the proximate lengths
    private int mGlobalShortInvCount;
    private int mGlobalLongDelDupInvCount;
    private double mArmLengthFactor;
    private final List<RegionAllocator> mDelRegionAllocators;
    private final List<RegionAllocator> mDupRegionAllocators;

    private boolean mLogVerbose;

    public static final int PRE_GENE_3P_DISTANCE = 10000;
    public static final int SHORT_INV_BUCKET = 100000;
    public static final int MIN_BUCKET_LENGTH = 10;
    public static final int LONG_DDI_BUCKET = 5000000;

    public static final double GENOME_BASE_COUNT = 3e9;
    public static final double MIN_FUSION_RATE = 1e-12;
    public static final String GENE_PAIR_DELIM = "_";

    public CohortExpFusions()
    {
        mChrGeneDataMap = Maps.newHashMap();
        mDelGenePairCounts = Maps.newHashMap();
        mDupGenePairCounts = Maps.newHashMap();
        mGlobalProximateCounts = Lists.newArrayList();
        mProximateBucketLengths = Lists.newArrayList();
        mMaxProximateBucketLength = 0;
        mGlobalShortInvCount = 0;
        mGlobalLongDelDupInvCount = 0;
        mArmLengthFactor = 0;
        mDelRegionAllocators = Lists.newArrayList();
        mDupRegionAllocators = Lists.newArrayList();
        mLogVerbose = false;
    }

    public final Map<String, List<GeneRangeData>> getChrGeneRangeDataMap() { return mChrGeneDataMap; }
    public Map<String,Map<Integer,Long>> getDelGenePairCounts() { return mDelGenePairCounts; }
    public Map<String,Map<Integer,Long>> getDupGenePairCounts() { return mDupGenePairCounts; }

    public double getArmLengthFactor() { return mArmLengthFactor; }

    public void initialiseLengths(final List<Integer> proximateBucketLengths, final List<String> restrictedChromosomes)
    {
        mProximateBucketLengths.addAll(proximateBucketLengths);

        mProximateBucketLengths.stream().forEach(x -> mGlobalProximateCounts.add(0));

        mMaxProximateBucketLength = !mProximateBucketLengths.isEmpty() ?
                mProximateBucketLengths.get(mProximateBucketLengths.size() - 1) : MIN_BUCKET_LENGTH;

        // sum up all arm lengths to adjust same-arm fusion rates
        int maxBucketLength = max(getMaxBucketLength(), 4000000);

        for(Map.Entry<Chromosome,Long> entry : SvUtilities.refGenomeLengths().lengths().entrySet())
        {
            final String chromosome = entry.getKey().toString();

            // if(!restrictedChromosomes.isEmpty() && !restrictedChromosomes.contains(chromosome))
            //    continue;

            int armLength = getChromosomalArmLength(chromosome, P_ARM);

            if(armLength > maxBucketLength)
                mArmLengthFactor += pow(armLength - maxBucketLength, 2);

            armLength = getChromosomalArmLength(chromosome, Q_ARM);

            if(armLength > maxBucketLength)
                mArmLengthFactor += armLength * (armLength - maxBucketLength);
        }

        FLC_LOGGER.debug("arm length factor: {}", mArmLengthFactor);

        // create region allocators for same-gene fusions
        int bucketLengths = mProximateBucketLengths.size() - 1;

        for(int i = 0; i < bucketLengths; ++i)
        {
            int blockSize = (int)(mProximateBucketLengths.get(i) / DEFAULT_BUCKET_REGION_RATIO);
            blockSize = max(blockSize, MIN_BLOCK_SIZE);
            mDelRegionAllocators.add(i, new RegionAllocator(blockSize));
            mDupRegionAllocators.add(i, new RegionAllocator(blockSize));
        }
    }

    public void setLogVerbose(boolean toggle) { mLogVerbose = toggle; }

    public void generateExpectedFusions(final EnsemblDataCache geneTransCache,
            final List<String> restrictedChromosomes, final List<String> restrictedGeneIds)
    {
        /* for each chromosome and arm:
            - convert each gene's transcripts into a set of phase regions
            - mark any overlapping phase regions between genes
            - determine same-gene fusions
            - determine proximate fusions
            - determine non-proximate fusions
            - finally determine global fusion rates
        */

        final Map<String, List<GeneData>> chrGeneDataMap = geneTransCache.getChrGeneDataMap();

        for(Map.Entry<String, List<GeneData>> entry : chrGeneDataMap.entrySet())
        {
            final String chromosome = entry.getKey();

            if(!restrictedChromosomes.isEmpty() && !restrictedChromosomes.contains(chromosome))
                continue;

            final List<GeneData> chrGenes = entry.getValue();

            if(chrGenes.isEmpty())
                continue;

            FLC_LOGGER.info("generating phase counts for chromosome({})", chromosome);

            List<GeneRangeData> chrGeneList = Lists.newArrayList();
            List<GeneRangeData> armGeneList = Lists.newArrayList();
            List<GeneRangeData> armGeneEndFirstList = Lists.newArrayList();

            ChromosomeArm currentArm = UNKNOWN;

            for(final GeneData geneData : chrGenes)
            {
                if(!restrictedGeneIds.isEmpty() && !restrictedGeneIds.contains(geneData.GeneId))
                    continue;

                GeneRangeData geneRangeData = new GeneRangeData(geneData);

                if(currentArm == UNKNOWN)
                {
                    currentArm = geneRangeData.Arm;
                }
                else if(currentArm != geneRangeData.Arm)
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
                    if(geneData.GeneEnd < armGeneEndFirstList.get(index).GeneData.GeneEnd)
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
                        FLC_LOGGER.debug("gene({}) {}", region.GeneId, region.toString());
                    }

                    // validity check
                    if(region.length() <= 0)
                    {
                        FLC_LOGGER.error("invalid region: gene({}) range({} -> {}) phase({})",
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

    private void processArmGenes(final String chromosome, final ChromosomeArm arm,
            List<GeneRangeData> armGeneList, List<GeneRangeData> armGeneEndFirstList)
    {
        FLC_LOGGER.info("chr({}) arm({}) processing {} genes", chromosome, arm, armGeneList.size());

        // mark all overlapping regions to avoid reallocation of base overlaps
        divideOverlappingRegions(chromosome, arm, armGeneList);

        for(int i = 0; i <= 1; ++i)
        {
            int strand = (i == 0) ? 1 : -1;

            List<GeneRangeData> geneDataList = strand == 1 ? armGeneList : armGeneEndFirstList;

            FLC_LOGGER.info("chr({}) arm({}) finding same-gene fusions", chromosome, arm);

            for(GeneRangeData geneData : geneDataList)
            {
                if(geneData.GeneData.Strand != strand)
                    continue;

                if(geneData.hasProteinCoding())
                {
                    generateSameGeneCounts(geneData);
                }
            }

            FLC_LOGGER.info("chr({}) arm({}) finding proximate fusions", chromosome, arm);

            generateProximateCounts(geneDataList, strand);

            FLC_LOGGER.info("chr({}) arm({}) finding non-proximate fusions", chromosome, arm);

            generateNonProximateCounts(geneDataList, strand);
        }

        generateNonProximateCounts(armGeneList, 0);
    }

    public void generatePhaseRegions(
            GeneRangeData geneRangeData, final List<TranscriptData> transDataList, final EnsemblDataCache geneTransCache)
    {
        // convert each transcript's exons into a set of phase regions spanning each intronic section
        List<GenePhaseRegion> phaseRegions = Lists.newArrayList();
        List<GenePhaseRegion> allTranscriptRegions = Lists.newArrayList();

        for(TranscriptData transcript : transDataList)
        {
            int transId = transcript.TransId;
            int precedingGeneSAPos = geneTransCache.findPrecedingGeneSpliceAcceptorPosition(transId);

            List<GenePhaseRegion> transcriptRegions = createPhaseRegionsFromTranscript(geneRangeData.GeneData, transcript, precedingGeneSAPos);

            transcriptRegions.forEach(x -> x.setTransId(transId));

            // add the new transcript's set of regions only if they aren't a very close overlap with any existing regions of the same phase
            for(GenePhaseRegion transRegion : transcriptRegions)
            {
                if(transRegion.hasPreGeneStatus() || !transRegion.hasPhasedType())
                    continue;

                if(!overlapsOtherRegions(transRegion, allTranscriptRegions, true, 0.90))
                    allTranscriptRegions.add(transRegion);
            }

            // consolidate regions where phases and pre-gene status overlap
            transcriptRegions.stream().forEach(x -> checkAddCombinedGenePhaseRegion(x, phaseRegions));
        }

        geneRangeData.setTranscriptPhaseRegions(allTranscriptRegions);

        mergePhaseRegions(phaseRegions);
        geneRangeData.setPhaseRegions(phaseRegions);
    }

    public static List<GenePhaseRegion> createPhaseRegionsFromTranscript(
            final GeneData geneData, final TranscriptData transcript, int precSpliceAcceptorPos)
    {
        // 5-prime rules - must be post-promotor
        // 3-prime rules: must be coding and > 1 exon, needs to find first splice acceptor and then uses its phasing
        // need to mark coding and non-coding regions

        List<GenePhaseRegion> transcriptRegions = Lists.newArrayList();

        // skip single-exon transcripts since without an intronic section their fusion likelihood is negligible
        if(transcript.exons().size() == 1)
            return transcriptRegions;

        boolean proteinCoding = transcript.BioType.equals(BIOTYPE_PROTEIN_CODING);

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

            if(geneData.Strand == 1 && nextExonData.Start > transcript.CodingEnd)
                break;
            else if(geneData.Strand == -1 && exonData.End < transcript.CodingStart)
                continue;

            // add an upstream gene region (only applicable for downstream fusion genes)
            // with a phase of -1 unless coding starts in the first exon or on the first base of the second exon
            if(precSpliceAcceptorPos > 0 && geneData.Strand == 1 && i == 0)
            {
                int regionPhase = (transcript.CodingStart < exonData.End) ? nextExonData.PhaseStart : -1;

                int preDistance = max(transcript.TransStart - precSpliceAcceptorPos, 0);
                int upstreamDistance = min(preDistance, PRE_GENE_3P_DISTANCE);
                int regionStart = transcript.TransStart - upstreamDistance;
                int regionEnd = exonData.Start - 1;

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
                int regionPhase = (transcript.CodingEnd > nextExonData.Start) ? exonData.PhaseStart : -1;

                int regionStart = nextExonData.End + 1;

                int preDistance = max(precSpliceAcceptorPos - transcript.TransEnd, 0);
                int upstreamDistance = min(preDistance, PRE_GENE_3P_DISTANCE);
                int regionEnd = transcript.TransEnd + upstreamDistance;

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
                        geneData.GeneId, exonData.Start, nextExonData.Start - 1, mapExonPhase(exonData.PhaseEnd));
            }
            else
            {
                phaseRegion = new GenePhaseRegion(
                        geneData.GeneId, exonData.End + 1, nextExonData.End, mapExonPhase(nextExonData.PhaseEnd));
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

        // clear any previous allocations from earlier genes
        mDelRegionAllocators.forEach(RegionAllocator::reset);
        mDupRegionAllocators.forEach(RegionAllocator::reset);

        for (int i = 0; i < transcriptRegions.size(); ++i)
        {
            GenePhaseRegion region1 = transcriptRegions.get(i);

            for (int j = i + 1; j < transcriptRegions.size(); ++j)
            {
                GenePhaseRegion region2 = transcriptRegions.get(j);

                if (region1.transId() != region2.transId())
                    break;

                testProximatePhaseRegions(geneData, geneData, region1, region2, true);
            }
        }
    }

    public int getMaxBucketLength() { return mMaxProximateBucketLength; }

    public void generateProximateCounts(final List<GeneRangeData> geneList, int strandMatch)
    {
        int rangeLimit = getMaxBucketLength();

        for(int lowerIndex = 0; lowerIndex < geneList.size(); ++lowerIndex)
        {
            GeneRangeData lowerGene = geneList.get(lowerIndex);

            if(lowerGene.GeneData.Strand != strandMatch)
                continue;

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
                        testProximatePhaseRegions(lowerGene, upperGene, lowerRegion, upperRegion, false);
                    }
                }
            }
        }
    }

    public void generateNonProximateCounts(final List<GeneRangeData> geneList, int strandMatch)
    {
        // test for candidate fusions between long DELs and DUPs, and then for INVs of various lengths
        int proximateLimit = getMaxBucketLength();

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
                int maxDistance = max(abs(gene1.GeneData.GeneStart - gene2.GeneData.GeneStart),
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

                        int regionOverlap = region1.length() * region2.length();

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

            FLC_LOGGER.info("calculating remote overlap counts for chr({})", chromosome);

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

    public static final int PERMITTED_REGION_OVERLAP = 5;

    public boolean testProximatePhaseRegions(
            GeneRangeData gene1, GeneRangeData gene2, GenePhaseRegion region1, GenePhaseRegion region2, boolean trackAllocations)
    {
        // ignore overlapping regions (allowing a small buffer) for now since it's not clear whether a DUP or DEL would be required
        if (haveOverlap(region1, region2, -PERMITTED_REGION_OVERLAP))
            return false;

        // skip tiny regions / overlaps
        if(region1.length() < PERMITTED_REGION_OVERLAP || region2.length() < PERMITTED_REGION_OVERLAP)
            return false;

        // test for DEL and DUP fusion independently
        boolean isForwardStrand = (gene1.GeneData.Strand == 1);
        boolean foundMatch = false;

        boolean region1IsLower = region1.start() < region2.start();
        GeneRangeData lowerGene = region1IsLower ? gene1 : gene2;
        GeneRangeData upperGene = region1IsLower ? gene2 : gene1;

        GenePhaseRegion lowerRegion = region1IsLower ? region1 :region2;
        GenePhaseRegion upperRegion = region1IsLower ? region2 :region1;

        for(int i = 0; i <= 1; ++i)
        {
            boolean isDel = (i == 0);

            boolean lowerGeneIsUpstream = (isDel == isForwardStrand);
            boolean upperGeneIsUpstream = !lowerGeneIsUpstream;

            if(upperGene.isStreamRestricted() && upperGeneIsUpstream != upperGene.isOnlyUpstreamPartner())
                continue;

            if(lowerGene.isStreamRestricted() && lowerGeneIsUpstream != lowerGene.isOnlyUpstreamPartner())
                continue;

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

            final Map<Integer,Long> bucketOverlapCounts = calcOverlapBucketAreas(
                    mProximateBucketLengths, trackAllocations ? (isDel ? mDelRegionAllocators : mDupRegionAllocators) : null,
                    lowerGene, upperGene, lowerRegion, upperRegion, isDel);

            if (!bucketOverlapCounts.isEmpty())
            {
                foundMatch = true;

                if (mLogVerbose)
                {
                    FLC_LOGGER.debug("gene({}: {}) and gene({}: {}) matched-region lower({} -> {}) upper({} -> {}) phase({}):",
                            lowerGene.GeneData.GeneId, lowerGene.GeneData.GeneName,
                            upperGene.GeneData.GeneId, upperGene.GeneData.GeneName,
                            lowerRegion.start(), lowerRegion.end(), upperRegion.start(), upperRegion.end(),
                            lowerRegion.getCombinedPhase());
                }

                for (Map.Entry<Integer,Long> entry : bucketOverlapCounts.entrySet())
                {
                    int bucketIndex = entry.getKey();
                    long overlap = entry.getValue();
                    addGeneFusionData(lowerGene, upperGene, overlap, isDel, bucketIndex);

                    if (mLogVerbose)
                    {
                        int bucketLen = mProximateBucketLengths.get(bucketIndex);
                        FLC_LOGGER.debug("matched-region bucketLength({}: {}) overlap({})", bucketIndex, bucketLen, overlap);
                    }
                }
            }
        }

        return foundMatch;
    }

    public static final int BUCKET_MIN = 0;
    public static final int BUCKET_MAX = 1;

    public int[] getBucketLengthMinMax(int bucketIndex)
    {
        if(bucketIndex >= mProximateBucketLengths.size())
            return new int[2];

        return new int[] {mProximateBucketLengths.get(bucketIndex), mProximateBucketLengths.get(bucketIndex + 1)};
    }

    public void addGeneFusionData(final GeneRangeData lowerGene, final GeneRangeData upperGene, long overlapCount, boolean isDel, int bucketIndex)
    {
        int[] bucketMinMax = getBucketLengthMinMax(bucketIndex);
        int bucketWidth = bucketMinMax[BUCKET_MAX] - bucketMinMax[BUCKET_MIN];
        double fusionRate = overlapCount / (bucketWidth * GENOME_BASE_COUNT);

        if(fusionRate < MIN_FUSION_RATE)
            return;

        if(mLogVerbose)
        {
            FLC_LOGGER.debug("gene pair({} & {}) adding {} overlap({}) for bucket length index({})",
                    lowerGene.GeneData.GeneName, upperGene.GeneData.GeneName,
                    isDel ? "DEL" : "DUP", overlapCount, bucketIndex);
        }

        Map<String,Map<Integer,Long>> genePairCounts = isDel ? mDelGenePairCounts : mDupGenePairCounts;

        final String genePair = lowerGene.GeneData.GeneId + GENE_PAIR_DELIM + upperGene.GeneData.GeneId;
        Map<Integer,Long> bucketOverlapCounts = genePairCounts.get(genePair);

        if(bucketOverlapCounts == null)
        {
            bucketOverlapCounts = Maps.newHashMap();
            genePairCounts.put(genePair, bucketOverlapCounts);
        }

        setBucketLengthData(bucketOverlapCounts, bucketIndex, overlapCount);
    }

    public final Map<String,GeneRangeData> generateGeneRangeData(final EnsemblDataCache geneTransCache, List<String> geneIds)
    {
        final Map<String, GeneRangeData> geneIdRangeDataMap = Maps.newHashMap();

        for(final String geneId : geneIds)
        {
            final GeneData geneData = geneTransCache.getGeneDataById(geneId);

            if(geneData == null)
            {
                FLC_LOGGER.error("geneId({}) not found", geneId);
                continue;
            }

            GeneRangeData geneRangeData = new GeneRangeData(geneData);

            // load from Ensembl transcript and exon data
            final List<TranscriptData> transDataList = geneTransCache.getTranscripts(geneData.GeneId);

            if (transDataList == null || transDataList.isEmpty())
                continue;

            generatePhaseRegions(geneRangeData, transDataList, geneTransCache);

            for(GenePhaseRegion region : geneRangeData.getPhaseRegions())
            {
                if(mLogVerbose)
                {
                    FLC_LOGGER.debug("gene({}) {}", region.GeneId, region.toString());
                }

                // validity check
                if(region.length() <= 0)
                {
                    FLC_LOGGER.error("invalid region: gene({}) range({} -> {}) phase({})",
                            region.GeneId, region.start(), region.end(), region.getCombinedPhase());
                    continue;
                }
            }

            geneIdRangeDataMap.put(geneId, geneRangeData);
        }

        return geneIdRangeDataMap;
    }

    public void logGlobalCounts()
    {
        FLC_LOGGER.info("GLOBAL_COUNTS: Type,GlobalCount,LengthMin,LengthMax");
        FLC_LOGGER.info("GLOBAL_COUNTS: LONG_DDI,{},{},{}", mGlobalLongDelDupInvCount, 5000000, 5000000);
        FLC_LOGGER.info("GLOBAL_COUNTS: INV,{},{},{}", mGlobalShortInvCount, 1000, 2000);

        // assign to correct bucket length
        for(int b = 0; b < mProximateBucketLengths.size() - 1; ++b)
        {
            int minLength = mProximateBucketLengths.get(b);
            int maxLength = mProximateBucketLengths.get(b + 1);

            FLC_LOGGER.info("GLOBAL_COUNTS: DEL,{},{},{}", mGlobalProximateCounts.get(b), minLength, maxLength);
            FLC_LOGGER.info("GLOBAL_COUNTS: DUP,{},{},{}", mGlobalProximateCounts.get(b), minLength, maxLength);
        }
    }


}
