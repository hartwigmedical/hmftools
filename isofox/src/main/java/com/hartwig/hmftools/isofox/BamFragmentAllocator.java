package com.hartwig.hmftools.isofox;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.fusion.KnownFusionType.KNOWN_PAIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.MULTI_MAP_QUALITY_THRESHOLD;
import static com.hartwig.hmftools.isofox.IsofoxConstants.SINGLE_MAP_QUALITY;
import static com.hartwig.hmftools.isofox.IsofoxFunction.NOVEL_LOCATIONS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.TRANSCRIPT_COUNTS;
import static com.hartwig.hmftools.isofox.IsofoxFunction.UNMAPPED_READS;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.DISCORDANT;
import static com.hartwig.hmftools.isofox.common.FragmentType.ALT;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.LOW_MAP_QUAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.isofox.common.FragmentType.UNSPLICED;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.common.ReadRecord.MAX_SC_BASE_MATCH;
import static com.hartwig.hmftools.isofox.common.ReadRecord.findOverlappingRegions;
import static com.hartwig.hmftools.isofox.common.ReadRecord.getUniqueValidRegion;
import static com.hartwig.hmftools.isofox.common.ReadRecord.markRegionBases;
import static com.hartwig.hmftools.isofox.common.ReadRecord.validTranscriptType;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_BOUNDARY;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.validExonMatch;
import static com.hartwig.hmftools.isofox.common.RnaUtils.deriveCommonRegions;
import static com.hartwig.hmftools.isofox.common.TransMatchType.OTHER_TRANS;
import static com.hartwig.hmftools.isofox.common.TransMatchType.SPLICE_JUNCTION;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;

import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.adjusts.GcRatioCounts.calcGcRatioFromReadRegions;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.isRealignedFragmentCandidate;
import static com.hartwig.hmftools.isofox.fusion.FusionUtils.setHasMultipleKnownSpliceGenes;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.fusion.KnownFusionData;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.DuplicateReadTracker;
import com.hartwig.hmftools.isofox.common.FragmentMatchType;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.FragmentType;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransMatchType;
import com.hartwig.hmftools.isofox.expression.CategoryCountsData;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.fusion.ChimericReadTracker;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder;
import com.hartwig.hmftools.isofox.novel.RetainedIntronFinder;
import com.hartwig.hmftools.isofox.novel.SpliceSiteCounter;
import com.hartwig.hmftools.isofox.results.ResultsWriter;
import com.hartwig.hmftools.isofox.unmapped.UnmappedReads;

import org.jetbrains.annotations.NotNull;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BamFragmentAllocator
{
    private final IsofoxConfig mConfig;
    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    // state relating to the current gene
    private GeneCollection mCurrentGenes;
    private final FragmentTracker mFragmentReads; // delay processing of read until both have been read

    private int mGeneReadCount;
    private int mTotalBamReadCount;
    private int mNextGeneCountLog;

    private final List<CategoryCountsData> mTransComboData;
    private final AltSpliceJunctionFinder mAltSpliceJunctionFinder;
    private final RetainedIntronFinder mRetainedIntronFinder;
    private final ChimericReadTracker mChimericReads;
    private final SpliceSiteCounter mSpliceSiteCounter;
    private final UnmappedReads mUnmappedReads;
    private final int[] mValidReadStartRegion;
    private final BaseDepth mBaseDepth;

    private final boolean mRunFusions;
    private final boolean mFusionsOnly;
    private final boolean mStatsOnly;

    private final BufferedWriter mReadDataWriter;
    private final GcRatioCounts mGcRatioCounts;
    private final GcRatioCounts mGeneGcRatioCounts;
    private int mEnrichedGeneFragments;
    private final DuplicateReadTracker mDuplicateTracker;
    private final List<String[]> mKnownPairGeneIds;
    private ChrBaseRegion mExcludedRegion;

    public BamFragmentAllocator(final IsofoxConfig config, final ResultsWriter resultsWriter)
    {
        mConfig = config;

        mCurrentGenes = null;
        mFragmentReads = new FragmentTracker();
        mTransComboData = Lists.newArrayList();

        mRunFusions = mConfig.Functions.contains(FUSIONS);
        mFusionsOnly = mConfig.runFusionsOnly();
        mStatsOnly = mConfig.runStatisticsOnly();

        mGeneReadCount = 0;
        mTotalBamReadCount = 0;
        mNextGeneCountLog = 0;
        mEnrichedGeneFragments = 0;
        mExcludedRegion = null;
        mValidReadStartRegion = new int[SE_PAIR];

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        // duplicates aren't counted towards fusions so can be ignored if only running fusions
        // reads with supplementary alignment data are only used for fusions
        boolean keepDuplicates = mConfig.runFunction(TRANSCRIPT_COUNTS);
        boolean keepSupplementaries = mRunFusions || mConfig.runFunction(NOVEL_LOCATIONS);
        boolean keepSecondaries = mConfig.ApplyMapQualityAdjust;
        int minMapQuality = keepSecondaries ? 0 : SINGLE_MAP_QUALITY;

        mBamSlicer = new BamSlicer(minMapQuality, keepDuplicates, keepSupplementaries, keepSecondaries);

        mDuplicateTracker = new DuplicateReadTracker(mConfig.MarkDuplicates);

        mReadDataWriter = resultsWriter.getReadDataWriter();
        mGcRatioCounts = mConfig.requireGcRatioCalcs() ? new GcRatioCounts() : null;
        mGeneGcRatioCounts = mConfig.requireGcRatioCalcs() ? new GcRatioCounts() : null;
        mBaseDepth = new BaseDepth();
        mChimericReads = new ChimericReadTracker(mConfig);
        mSpliceSiteCounter = new SpliceSiteCounter(resultsWriter.getSpliceSiteWriter());

        mKnownPairGeneIds = Lists.newArrayList();

        mAltSpliceJunctionFinder = new AltSpliceJunctionFinder(mConfig, resultsWriter.getAltSpliceJunctionWriter());
        mRetainedIntronFinder = new RetainedIntronFinder(mConfig, resultsWriter.getRetainedIntronWriter());
        mUnmappedReads = mConfig.runFunction(UNMAPPED_READS) ? new UnmappedReads(mConfig, resultsWriter.getUnmappedReadsWriter()) : null;
    }

    public int totalReadCount() { return mTotalBamReadCount; }

    public final GcRatioCounts getGcRatioCounts() { return mGcRatioCounts; }
    public final GcRatioCounts getGeneGcRatioCounts() { return mGeneGcRatioCounts; }
    public BaseDepth getBaseDepth() { return mBaseDepth; }
    public final ChimericReadTracker getChimericReadTracker() { return mChimericReads; }
    public final SpliceSiteCounter getSpliceSiteCounter() { return mSpliceSiteCounter; }

    private static int GENE_LOG_COUNT = 5000000;

    public void clearCache()
    {
        mFragmentReads.clear();
        mTransComboData.clear();
        mChimericReads.clear();

        if(mGeneGcRatioCounts != null)
            mGeneGcRatioCounts.clearCounts();

        mAltSpliceJunctionFinder.setGeneData(null);
        mRetainedIntronFinder.setGeneData(null);
        mSpliceSiteCounter.clear();

        mDuplicateTracker.clear();
        mCurrentGenes = null;
        mExcludedRegion = null;
    }

    private static final int NON_GENIC_BASE_DEPTH_WIDTH = 250000;

    public void produceBamCounts(final GeneCollection geneCollection, final ChrBaseRegion geneRegion)
    {
        clearCache();

        mCurrentGenes = geneCollection;

        mGeneReadCount = 0;
        mNextGeneCountLog = GENE_LOG_COUNT;
        mEnrichedGeneFragments = 0;

        // and width around the base depth region to pick up junctions outside the gene
        int[] baseDepthRange = new int[SE_PAIR];
        baseDepthRange[SE_START] = max(geneRegion.start(), geneCollection.regionBounds()[SE_START] - NON_GENIC_BASE_DEPTH_WIDTH);
        baseDepthRange[SE_END] = min(geneRegion.end(), geneCollection.regionBounds()[SE_END] + NON_GENIC_BASE_DEPTH_WIDTH);
        mBaseDepth.initialise(baseDepthRange);

        mChimericReads.initialise(mCurrentGenes);

        if(mAltSpliceJunctionFinder.enabled())
            mAltSpliceJunctionFinder.setGeneData(mCurrentGenes);

        if(mRetainedIntronFinder.enabled())
            mRetainedIntronFinder.setGeneData(mCurrentGenes);

        mValidReadStartRegion[SE_START] = geneRegion.start();
        mValidReadStartRegion[SE_END] = geneRegion.end();

        if(mConfig.ExcludedRegion != null && geneRegion.overlaps(mConfig.ExcludedRegion))
        {
            // special handling to avoid any specified enriched region (in this case LINC00486's poly-G sequence)
            mExcludedRegion = mConfig.ExcludedRegion;
            final ChrBaseRegion preRegion = new ChrBaseRegion(geneRegion.Chromosome, geneRegion.start(), mExcludedRegion.start() - 100);
            final ChrBaseRegion postRegion = new ChrBaseRegion(geneRegion.Chromosome, mExcludedRegion.end() + 100, geneRegion.end());
            mBamSlicer.slice(mSamReader, Lists.newArrayList(preRegion), this::processSamRecord);
            mBamSlicer.slice(mSamReader, Lists.newArrayList(postRegion), this::processSamRecord);
        }
        else
        {
            mBamSlicer.slice(mSamReader, Lists.newArrayList(geneRegion), this::processSamRecord);
        }

        if(mEnrichedGeneFragments > 0)
            processEnrichedGeneFragments();

        if(mChimericReads.enabled())
        {
            mChimericReads.postProcessChimericReads(mBaseDepth, mFragmentReads);
            processChimericNovelJunctions();
        }

        ISF_LOGGER.debug("genes({}) bamReadCount({}) depth(bases={} perc={} max={})",
                mCurrentGenes.geneNames(), mGeneReadCount, mBaseDepth.basesWithDepth(),
                String.format("%.3f", mBaseDepth.basesWithDepthPerc()), mBaseDepth.maxDepth());
    }

    private void processChimericNovelJunctions()
    {
        if(!mAltSpliceJunctionFinder.enabled() || mChimericReads.getLocalChimericReads().isEmpty())
            return;

        final List<Integer> invalidTrans = Lists.newArrayList();

        for(final List<ReadRecord> reads : mChimericReads.getLocalChimericReads())
        {
            ReadRecord read1 = null;
            ReadRecord read2 = null;

            if(reads.size() == 2)
            {
                read1 = reads.get(0);
                read2 = reads.get(1);
            }
            else if(reads.size() == 3)
            {
                for(ReadRecord read : reads)
                {
                    if(read.hasSuppAlignment())
                    {
                        if(read1 == null)
                        {
                            read1 = read;
                        }
                        else
                        {
                            read2 = read;
                            break;
                        }
                    }
                }
            }

            if(read1 == null || read2 == null)
                continue;

            int readPosMin = min(read1.PosStart, read2.PosStart);
            int readPosMax = max(read1.PosEnd, read2.PosEnd);

            final List<GeneReadData> overlapGenes = mCurrentGenes.findGenesCoveringRange(readPosMin, readPosMax, false);
            mAltSpliceJunctionFinder.evaluateFragmentReads(overlapGenes, read1, read2, invalidTrans);
        }
    }

    public void annotateNovelLocations()
    {
        recordNovelLocationReadDepth();
        mAltSpliceJunctionFinder.prioritiseGenes();
        mAltSpliceJunctionFinder.writeAltSpliceJunctions();
        mRetainedIntronFinder.writeRetainedIntrons();
    }

    private void processSamRecord(@NotNull final SAMRecord record)
    {
        // to avoid double-processing of reads overlapping 2 (or more) gene collections, only process them if they start in this
        // gene collection or its preceding non-genic region
        if(!positionWithin(record.getStart(), mValidReadStartRegion[SE_START], mValidReadStartRegion[SE_END]))
            return;

        if(excludeRegion(record))
            return;

        if(mCurrentGenes.inEnrichedRegion(record.getStart(), record.getEnd()))
        {
            processEnrichedRegionRead(record);
            return;
        }

        if(mDuplicateTracker.checkDuplicates(record))
        {
            if(mConfig.DropDuplicates)
            {
                if(record.getFirstOfPairFlag())
                    mCurrentGenes.addCount(DUPLICATE, 1);

                return;
            }
        }

        ++mTotalBamReadCount;
        ++mGeneReadCount;

        processRead(ReadRecord.from(record));
    }

    // private static final String LOG_READ_ID = "";
    private static final String LOG_READ_ID = "";

    private void processRead(ReadRecord read)
    {
        if(read.Id.equals(LOG_READ_ID))
        {
            ISF_LOGGER.debug("specific read by ID: {}", read.toString());
        }

        if(mUnmappedReads != null)
            mUnmappedReads.processReadRecord(read, mCurrentGenes);

        // for each record find all exons with an overlap
        // skip records if either end isn't in one of the exons for this gene

        if(mGeneReadCount >= mNextGeneCountLog)
        {
            mNextGeneCountLog += GENE_LOG_COUNT;
            ISF_LOGGER.info("chr({}) genes({}) bamRecordCount({})", mCurrentGenes.chromosome(), mCurrentGenes.geneNames(), mGeneReadCount);
        }

        if(mConfig.GeneReadLimit > 0 && mGeneReadCount >= mConfig.GeneReadLimit)
        {
            mBamSlicer.haltProcessing();
            ISF_LOGGER.warn("chr({}) genes({}) readCount({}) exceeds max read count",
                    mCurrentGenes.chromosome(), mCurrentGenes.geneNames(), mGeneReadCount);
            return;
        }

        final List<RegionReadData> overlappingRegions = findOverlappingRegions(mCurrentGenes.getExonRegions(), read);

        if (!overlappingRegions.isEmpty())
        {
            read.processOverlappingRegions(overlappingRegions);
        }

        mCurrentGenes.setReadGeneCollections(read, mValidReadStartRegion);

        checkFragmentRead(read);
    }

    private boolean checkFragmentRead(ReadRecord read)
    {
        // check if the 2 reads from a fragment exist and if so handle them a pair, returning true
        ReadRecord otherRead = mFragmentReads.checkRead(read);

        if(otherRead != null)
        {
            processFragmentReads(read, otherRead);
            return true;
        }

        return false;
    }

    private void processFragmentReads(final ReadRecord read1, final ReadRecord read2)
    {
        /* process the pair of reads from a fragment:
            - fully outside the gene (due to the buffer used, ignore
            - read through a gene ie start or end outside
            - purely intronic
            - chimeric an inversion or translocation
            - supporting 1 or more transcripts
                - both reads fully with an exon - if exon has only 1 transcript then consider unambiguous
                - both reads within 2 exons (including spanning intermediary ones) and/or either exon at the boundary
            - not supporting any transcript - eg alternative splice sites or unspliced reads
        */

        boolean isDuplicate = read1.isDuplicate() || read2.isDuplicate();
        int minMapQuality = min(read1.mapQuality(), read2.mapQuality());
        boolean isMultiMapped = minMapQuality <= MULTI_MAP_QUALITY_THRESHOLD;
        boolean isChimeric = read1.isChimeric() || read2.isChimeric() || !read1.withinGeneCollection() || !read2.withinGeneCollection();

        if(!isChimeric && !isDuplicate && !isMultiMapped && mChimericReads.enabled() && (read1.containsSplit() || read2.containsSplit()))
        {
            isChimeric = setHasMultipleKnownSpliceGenes(Lists.newArrayList(read1, read2), mKnownPairGeneIds);
        }

        if(mStatsOnly)
        {
            if(isDuplicate)
                mCurrentGenes.addCount(DUPLICATE, 1);
            else if(isChimeric)
                mCurrentGenes.addCount(CHIMERIC, 1);
            else if(read1.getMappedRegions().isEmpty() && read2.getMappedRegions().isEmpty())
                mCurrentGenes.addCount(UNSPLICED, 1);
            else
                mCurrentGenes.addCount(TRANS_SUPPORTING, 1);

            mCurrentGenes.addCount(TOTAL, 1);
            return;
        }

        final List<int[]> commonMappings = deriveCommonRegions(read1.getMappedRegionCoords(), read2.getMappedRegionCoords());

        if(isDuplicate)
        {
            mCurrentGenes.addCount(DUPLICATE, 1);
        }
        else
        {
            mBaseDepth.processRead(commonMappings);
        }

        // if either read is chimeric (including one outside the genic region) then handle them both as such
        // some of these may be re-processed as alternative SJ candidates if they are within a single gene
        if(isChimeric)
        {
            if(!isMultiMapped)
                processChimericReadPair(read1, read2);

            return;
        }

        if(mRunFusions)
        {
            // reads with sufficient soft-clipping and not mapped to an adjacent region are candidates for fusion re-alignment
            if(isRealignedFragmentCandidate(read1) || isRealignedFragmentCandidate(read2))
            {
                mChimericReads.addRealignmentCandidates(read1, read2);
            }

            if (mFusionsOnly)
                return;
        }

        int readPosMin = min(read1.PosStart, read2.PosStart);
        int readPosMax = max(read1.PosEnd, read2.PosEnd);

        final List<GeneReadData> overlapGenes = mCurrentGenes.findGenesCoveringRange(readPosMin, readPosMax, true);
        mCurrentGenes.addCount(TOTAL, 1);

        if(read1.getMappedRegions().isEmpty() && read2.getMappedRegions().isEmpty())
        {
            // fully intronic read in every transcript and gene
            processIntronicReads(overlapGenes, read1, read2);
            return;
        }

        final Map<Integer,TransMatchType> firstReadTransTypes = read1.getTranscriptClassifications();
        final Map<Integer,TransMatchType> secondReadTransTypes = read2.getTranscriptClassifications();

        // first find valid transcripts in both reads
        final List<Integer> validTranscripts = Lists.newArrayList();
        final Set<Integer> invalidTranscripts = Sets.newHashSet();

        final List<RegionReadData> validRegions = getUniqueValidRegion(read1, read2);

        if(mConfig.RunValidations)
        {
            for(RegionReadData region : validRegions)
            {
                if(validRegions.stream().filter(x -> x == region).count() > 1)
                {
                    ISF_LOGGER.error("repeated exon region({})", region);
                }
            }
        }

        // track splice site info
        if(mConfig.WriteSpliceSiteData)
        {
            mSpliceSiteCounter.registerSpliceSiteSupport(read1.getMappedRegionCoords(), read2.getMappedRegionCoords(), mCurrentGenes.getExonRegions());
        }

        for(Map.Entry<Integer,TransMatchType> entry : firstReadTransTypes.entrySet())
        {
            int transId = entry.getKey();

            if(validTranscriptType(entry.getValue()) && secondReadTransTypes.containsKey(transId) && validTranscriptType(secondReadTransTypes.get(transId)))
            {
                int calcFragmentLength = calcFragmentLength(transId, read1, read2);
                boolean validFragmentLength = calcFragmentLength > 0 && calcFragmentLength <= mConfig.MaxFragmentLength;

                if(validFragmentLength)
                {
                    validTranscripts.add(transId);
                }
                else
                {
                    invalidTranscripts.add(transId);
                }
            }
            else
            {
                invalidTranscripts.add(transId);
            }
        }

        for(Integer transId : secondReadTransTypes.keySet())
        {
            if(!validTranscripts.contains(transId))
                invalidTranscripts.add(transId);
        }

        FragmentType fragmentType = UNSPLICED;

        // now mark all other transcripts which aren't valid either due to the read pair
        if(validTranscripts.isEmpty())
        {
            // no valid transcripts but record against the gene further information about these reads
            boolean checkRetainedIntrons = false;

            if(read1.containsSplit() || read2.containsSplit())
            {
                fragmentType = ALT;

                if(mAltSpliceJunctionFinder.enabled())
                {
                    mAltSpliceJunctionFinder.evaluateFragmentReads(
                            overlapGenes, read1, read2, invalidTranscripts.stream().collect(Collectors.toList()));
                }

                checkRetainedIntrons = true;
            }
            else
            {
                // look for alternative splicing from long reads involving more than one region and not spanning into an intron
                for(int transId : invalidTranscripts)
                {
                    List<RegionReadData> regions = read1.getMappedRegions().entrySet().stream()
                            .filter(x -> x.getKey().hasTransId(transId))
                            .filter(x -> x.getValue() != EXON_INTRON)
                            .map(x -> x.getKey()).collect(Collectors.toList());;

                    final List<RegionReadData> regions2 = read2.getMappedRegions().entrySet().stream()
                            .filter(x -> x.getKey().hasTransId(transId))
                            .filter(x -> x.getValue() != EXON_INTRON)
                            .map(x -> x.getKey()).collect(Collectors.toList());

                    for(RegionReadData region : regions2)
                    {
                        if (!regions.contains(region))
                            regions.add(region);
                    }

                    if(regions.size() > 1)
                    {
                        fragmentType = ALT;
                        break;
                    }
                }

                checkRetainedIntrons = true;
            }

            if(checkRetainedIntrons && mRetainedIntronFinder.enabled())
                mRetainedIntronFinder.evaluateFragmentReads(read1, read2);

            if(fragmentType == UNSPLICED)
            {
                List<String> unsplicedGeneIds = overlapGenes.stream().map(x -> x.GeneData.GeneId).collect(Collectors.toList());

                if(!unsplicedGeneIds.isEmpty())
                {
                    CategoryCountsData catCounts = getCategoryCountsData(validTranscripts, unsplicedGeneIds);
                    addGcCounts(catCounts, commonMappings, minMapQuality);
                }
            }
        }
        else
        {
            // record valid read info against each region now that it is known
            fragmentType = TRANS_SUPPORTING;

            // first mark any invalid trans as 'other' meaning it doesn't require any further classification since a valid trans exists
            firstReadTransTypes.entrySet().stream()
                    .filter(x -> validTranscriptType(x.getValue()))
                    .filter(x -> !validTranscripts.contains(x.getKey()))
                    .forEach(x -> x.setValue(OTHER_TRANS));

            secondReadTransTypes.entrySet().stream()
                    .filter(x -> validTranscriptType(x.getValue()))
                    .filter(x -> !validTranscripts.contains(x.getKey()))
                    .forEach(x -> x.setValue(OTHER_TRANS));

            if(mConfig.RunValidations)
            {
                for(int[] readRegion : commonMappings)
                {
                    if(commonMappings.stream().filter(x -> x[SE_START] == readRegion[SE_START] && x[SE_END] == readRegion[SE_END]).count() > 1)
                    {
                        ISF_LOGGER.error("repeated read region({} -> {})", readRegion[SE_START], readRegion[SE_END]);
                    }
                }
            }

            validRegions.forEach(x -> markRegionBases(commonMappings, x));

            // now set counts for each valid transcript
            boolean isUniqueTrans = validTranscripts.size() == 1;

            FragmentMatchType comboTransMatchType = FragmentMatchType.SHORT;

            for (int transId : validTranscripts)
            {
                int regionCount = (int)validRegions.stream().filter(x -> x.hasTransId(transId)).count();

                FragmentMatchType transMatchType;

                if(read1.getTranscriptClassification(transId) == SPLICE_JUNCTION || read2.getTranscriptClassification(transId) == SPLICE_JUNCTION)
                {
                    transMatchType = FragmentMatchType.SPLICED;
                    comboTransMatchType = FragmentMatchType.SPLICED;
                }
                else if(regionCount > 1)
                {
                    // the read pair span at least 2 exons within the same transcript
                    transMatchType = FragmentMatchType.LONG;

                    if(comboTransMatchType != FragmentMatchType.SPLICED)
                        comboTransMatchType = FragmentMatchType.LONG;
                }
                else
                {
                    transMatchType = FragmentMatchType.SHORT;
                }

                mCurrentGenes.addTranscriptReadMatch(transId, isUniqueTrans, transMatchType);

                // separately record discordant reads spanning 2+ exons
                if(!read1.containsSplit() && read2.containsSplit())
                {
                    final RegionReadData region1 =
                            read1.getMappedRegions().keySet().stream().filter(x -> x.hasTransId(transId)).findFirst().orElse(null);
                    final RegionReadData region2 =
                            read2.getMappedRegions().keySet().stream().filter(x -> x.hasTransId(transId)).findFirst().orElse(null);

                    if(region1 != region2 && region1.getExonRank(transId) != region2.getExonRank(transId))
                    {
                        mCurrentGenes.addTranscriptReadMatch(transId, DISCORDANT);
                    }
                }

                // keep track of which regions have been allocated from this fragment as a whole, so not counting each read separately
                processValidTranscript(transId, Lists.newArrayList(read1, read2), isUniqueTrans);
            }

            List<String> unsplicedGeneIds = comboTransMatchType == FragmentMatchType.SHORT ?
                    overlapGenes.stream().map(x -> x.GeneData.GeneId).collect(Collectors.toList()) : Lists.newArrayList();

            CategoryCountsData catCounts = getCategoryCountsData(validTranscripts, unsplicedGeneIds);
            addGcCounts(catCounts, commonMappings, minMapQuality);
        }

        mCurrentGenes.addCount(fragmentType, 1);

        if(mConfig.WriteReadData && mReadDataWriter != null)
        {
            for(final GeneReadData geneReadData : overlapGenes)
            {
                writeReadData(mReadDataWriter, geneReadData, 0, read1, read2, fragmentType, validTranscripts.size());
                writeReadData(mReadDataWriter, geneReadData, 1, read2, read1, fragmentType, validTranscripts.size());
            }
        }
    }

    private int calcFragmentLength(int transId, final ReadRecord read1, final ReadRecord read2)
    {
        final TranscriptData transData = mCurrentGenes.getTranscripts().stream().filter(x -> x.TransId == transId).findFirst().orElse(null);
        if(transData == null)
            return -1;

        return calcFragmentLength(transData, read1, read2);
    }

    public static int calcFragmentLength(final TranscriptData transData, final ReadRecord read1, final ReadRecord read2)
    {
        int minReadPos = min(read1.PosStart, read2.PosStart);
        int maxReadPos = max(read1.PosEnd, read2.PosEnd);
        return calcFragmentLength(transData, minReadPos, maxReadPos);
    }

    public static int calcFragmentLength(final TranscriptData transData, final int minReadPos, final int maxReadPos)
    {
        // calculate fragment length within this transcript assuming it has been spliced
        int transcriptBases = 0;
        boolean startFound = false;

        for(final ExonData exon : transData.exons())
        {
            if(!startFound)
            {
                if(minReadPos < exon.Start - MAX_SC_BASE_MATCH)
                    break;

                if(minReadPos > exon.End)
                    continue;

                if(maxReadPos <= exon.End)
                {
                    // within same exon
                    return maxReadPos - minReadPos + 1;
                }

                startFound = true;
                transcriptBases = exon.End - max(exon.Start, minReadPos) + 1;
            }
            else
            {
                if(maxReadPos > exon.End)
                {
                    transcriptBases += exon.baseLength();
                }
                else if(maxReadPos < exon.Start)
                {
                    break;
                }
                else
                {
                    transcriptBases += maxReadPos - exon.Start + 1;
                    break;
                }
            }
        }

        return transcriptBases;
    }

    private boolean excludeRegion(final SAMRecord record)
    {
        if(mExcludedRegion == null)
            return false;

        return (mExcludedRegion.containsPosition(record.getStart()) || mExcludedRegion.containsPosition(record.getEnd())
                || mExcludedRegion.containsPosition(record.getMateAlignmentStart()));
    }

    private void processEnrichedRegionRead(final SAMRecord record)
    {
        ++mTotalBamReadCount;
        ++mGeneReadCount;

        if(mGeneReadCount >= mNextGeneCountLog)
        {
            mNextGeneCountLog += GENE_LOG_COUNT;
            ISF_LOGGER.info("chr({}) genes({}) bamRecordCount({})", mCurrentGenes.chromosome(), mCurrentGenes.geneNames(), mGeneReadCount);
        }

        if(mConfig.GeneReadLimit > 0 && mGeneReadCount >= mConfig.GeneReadLimit)
        {
            mBamSlicer.haltProcessing();
            ISF_LOGGER.warn("chr({}) genes({}) readCount({}) exceeds max read count",
                    mCurrentGenes.chromosome(), mCurrentGenes.geneNames(), mGeneReadCount);
            return;
        }

        if(record.getFirstOfPairFlag()) // only count once per fragment
            return;

        if(!mConfig.runFunction(TRANSCRIPT_COUNTS))
            return;

        if(record.getDuplicateReadFlag())
            mCurrentGenes.addCount(DUPLICATE, 1);

        // no further classification of fragment is performed - ie they are considered supporting
        ++mEnrichedGeneFragments;
    }

    private void processEnrichedGeneFragments()
    {
        // add to overall counts - since these are within a single exon, consider them supporting the transcript + unspliced
        mCurrentGenes.addCount(TOTAL, mEnrichedGeneFragments);
        mCurrentGenes.addCount(TRANS_SUPPORTING, mEnrichedGeneFragments);

        // add to category counts
        final int[] enrichedRegion = mCurrentGenes.getEnrichedRegion();
        final List<String> unsplicedGeneIds = mCurrentGenes.findGenesCoveringRange(enrichedRegion[SE_START], enrichedRegion[SE_END], true)
                .stream().map(x -> x.GeneData.GeneId).collect(Collectors.toList());

        final List<Integer> transIds = mCurrentGenes.getEnrichedTranscripts().stream().map(x -> new Integer(x.TransId)).collect(Collectors.toList());
        CategoryCountsData catCounts = getCategoryCountsData(transIds, unsplicedGeneIds);

        // compute and cache GC data
        double gcRatio = calcGcRatioFromReadRegions(mConfig.RefGenome, mCurrentGenes.chromosome(), Lists.newArrayList(mCurrentGenes.getEnrichedRegion()));

        int[] gcRatioIndices = { -1, -1 };
        double[] gcRatioCounts = { 0, 0 };

        if(mGcRatioCounts != null)
        {
            mGcRatioCounts.determineRatioData(gcRatio, gcRatioIndices, gcRatioCounts);
            gcRatioCounts[0] *= mEnrichedGeneFragments;
            gcRatioCounts[1] *= mEnrichedGeneFragments;
        }

        addGcCounts(catCounts, gcRatioIndices, gcRatioCounts, mEnrichedGeneFragments);
    }

    public List<CategoryCountsData> getTransComboData() { return mTransComboData; }

    private CategoryCountsData getCategoryCountsData(final List<Integer> transcripts, final List<String> geneIds)
    {
        CategoryCountsData transComboCounts = mTransComboData.stream()
                .filter(x -> x.matches(transcripts, geneIds)).findFirst().orElse(null);

        if(transComboCounts == null)
        {
            transComboCounts = new CategoryCountsData(transcripts, geneIds);

            if(mGcRatioCounts != null && mConfig.ApplyGcBiasAdjust)
                transComboCounts.initialiseGcRatioCounts(mGcRatioCounts.getCounts().length);

            mTransComboData.add(transComboCounts);
        }

        return transComboCounts;
    }

    private void processValidTranscript(int transId, final List<ReadRecord> reads, boolean isUniqueTrans)
    {
        final List<RegionReadData> processedRegions = Lists.newArrayList();

        for(ReadRecord read : reads)
        {
            List<RegionReadData> regions = read.getMappedRegions().entrySet().stream()
                    .filter(x -> x.getKey().hasTransId(transId))
                    .filter(x -> validExonMatch(x.getValue()))
                    .map(x -> x.getKey()).collect(Collectors.toList());

            for(RegionReadData region : regions)
            {
                if(!processedRegions.contains(region))
                {
                    // register a read against this valid transcript region
                    region.addTranscriptReadMatch(transId, isUniqueTrans);
                }
            }

            // any adjacent reads can record a splice junction count
            if(regions.size() > 1 && read.getTranscriptClassification(transId) == SPLICE_JUNCTION)
            {
                for(int r1 = 0; r1 < regions.size() - 1; ++r1)
                {
                    RegionReadData region1 = regions.get(r1);

                    for(int r2 = r1 + 1; r2 < regions.size(); ++r2)
                    {
                        RegionReadData region2 = regions.get(r2);

                        if(processedRegions.contains(region1) && processedRegions.contains(region2))
                            continue;

                        if(region1.getPostRegions().contains(region2))
                        {
                            region1.addTranscriptJunctionMatch(transId, SE_END, isUniqueTrans);
                            region2.addTranscriptJunctionMatch(transId, SE_START, isUniqueTrans);
                        }
                        else if(region1.getPreRegions().contains(region2))
                        {
                            region1.addTranscriptJunctionMatch(transId, SE_START, isUniqueTrans);
                            region2.addTranscriptJunctionMatch(transId, SE_END, isUniqueTrans);
                        }
                    }
                }
            }

            regions.forEach(x -> processedRegions.add(x));
        }
    }

    private void processIntronicReads(final List<GeneReadData> genes, final ReadRecord read1, final ReadRecord read2)
    {
        if(read1.Cigar.containsOperator(CigarOperator.N) || read2.Cigar.containsOperator(CigarOperator.N))
        {
            mCurrentGenes.addCount(ALT, 1);

            if(mAltSpliceJunctionFinder.enabled())
                mAltSpliceJunctionFinder.evaluateFragmentReads(genes, read1, read2, Lists.newArrayList());

            return;
        }

        List<String> unsplicedGeneIds = genes.stream().map(x -> x.GeneData.GeneId).collect(Collectors.toList());

        if(!unsplicedGeneIds.isEmpty())
        {
            CategoryCountsData catCounts = getCategoryCountsData(Lists.newArrayList(), unsplicedGeneIds);

            List<int[]> readRegions = deriveCommonRegions(read1.getMappedRegionCoords(), read2.getMappedRegionCoords());
            addGcCounts(catCounts, readRegions, min(read1.mapQuality(), read2.mapQuality()));
        }

        mCurrentGenes.addCount(UNSPLICED, 1);
    }

    private void addGcCounts(final CategoryCountsData catCounts, final List<int[]> readRegions, int minMapQuality)
    {
        int[] gcRatioIndices = { -1, -1 };
        double[] gcRatioCounts = { 0, 0 };

        if (mGcRatioCounts != null)
        {
            double gcRatio = calcGcRatioFromReadRegions(mConfig.RefGenome, mCurrentGenes.chromosome(), readRegions);
            mGcRatioCounts.determineRatioData(gcRatio, gcRatioIndices, gcRatioCounts);
        }

        double fragmentCount = 1;

        if(minMapQuality <= MULTI_MAP_QUALITY_THRESHOLD && mConfig.ApplyMapQualityAdjust)
        {
            if(minMapQuality == 3)
                fragmentCount = 0.5;
            else if(minMapQuality == 2)
                fragmentCount = 0.33;
            else if(minMapQuality == 1)
                fragmentCount = 0.2;
            else
                fragmentCount = 0.1;

            mCurrentGenes.addCount(LOW_MAP_QUAL, 1);
        }

        addGcCounts(catCounts, gcRatioIndices, gcRatioCounts, fragmentCount);
    }

    private void addGcCounts(final CategoryCountsData catCounts, final int[] gcRatioIndices, double[] gcRatioCounts, double count)
    {
        if(mGcRatioCounts != null)
        {
            for(int i = 0; i < gcRatioIndices.length; ++i)
            {
                if (gcRatioIndices[i] >= 0)
                {
                    mGcRatioCounts.addGcRatioCount(gcRatioIndices[i], gcRatioCounts[i]);
                    mGeneGcRatioCounts.addGcRatioCount(gcRatioIndices[i], gcRatioCounts[i]);
                }
            }

            if(mConfig.ApplyGcBiasAdjust)
                catCounts.addGcRatioCounts(count, gcRatioIndices, gcRatioCounts);
            else
                catCounts.addCounts(count);
        }
        else
        {
            catCounts.addCounts(count);
        }
    }

    private void recordNovelLocationReadDepth()
    {
        if(mAltSpliceJunctionFinder.getAltSpliceJunctions().isEmpty() && mRetainedIntronFinder.getRetainedIntrons().isEmpty())
            return;

        mAltSpliceJunctionFinder.setPositionDepth(mBaseDepth);
        mRetainedIntronFinder.setPositionDepth(mBaseDepth);
    }

    private void processChimericReadPair(final ReadRecord read1, final ReadRecord read2)
    {
        if(mChimericReads.enabled())
        {
            mChimericReads.addChimericReadPair(read1, read2);
        }
        else
        {
            // avoid double-counting fragment reads
            mCurrentGenes.addCount(TOTAL, 1);
            mCurrentGenes.addCount(CHIMERIC, 1);
        }
    }

    public void registerKnownFusionPairs(final EnsemblDataCache geneTransCache)
    {
        for(final KnownFusionData knownPair : mConfig.Fusions.KnownFusions.getDataByType(KNOWN_PAIR))
        {
            final GeneData upGene = geneTransCache.getGeneDataByName(knownPair.FiveGene);
            final GeneData downGene = geneTransCache.getGeneDataByName(knownPair.ThreeGene);
            if(upGene != null && downGene != null)
                mKnownPairGeneIds.add(new String[] { upGene.GeneId, downGene.GeneId });
        }

        mChimericReads.addKnownPairGeneIds(mKnownPairGeneIds);
    }

    public static BufferedWriter createReadDataWriter(final IsofoxConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile("read_data.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneId,GeneName,ReadIndex,ReadId,Chromosome,PosStart,PosEnd,Cigar");
            writer.write(",InsertSize,FragLength,MateChr,MatePosStart,FirstInPair,ReadReversed,SuppData");
            writer.write(",GeneClass,TransId,TransClass,ValidTrans,ExonRank,ExonStart,ExonEnd,RegionClass,ScRegionsStart,SvRegionsEnd");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to create read data writer: {}", e.toString());
            return null;
        }
    }

    private synchronized static void writeReadData(
            final BufferedWriter writer, final GeneReadData geneReadData, int readIndex, final ReadRecord read, final ReadRecord otherRead,
            FragmentType geneReadType, int validTranscripts)
    {
        if(read.getTranscriptClassifications().isEmpty())
            return;

        try
        {
            for(Map.Entry<Integer,TransMatchType> entry : read.getTranscriptClassifications().entrySet())
            {
                int transId = entry.getKey();
                TransMatchType transType = entry.getValue();

                for(Map.Entry<RegionReadData, RegionMatchType> rEntry : read.getMappedRegions().entrySet())
                {
                    RegionReadData region = rEntry.getKey();
                    RegionMatchType matchType = rEntry.getValue();

                    if(!region.hasTransId(transId))
                        continue;

                    writer.write(String.format("%s,%s,%d,%s",
                            geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName, readIndex, read.Id));

                    int calcFragmentLength = read.fragmentInsertSize();

                    if(validTranscriptType(read.getTranscriptClassification(transId)))
                    {
                        final TranscriptData transData =
                                geneReadData.getTranscripts().stream().filter(x -> x.TransId == transId).findFirst().orElse(null);
                        calcFragmentLength = transData != null ? calcFragmentLength(transData, read, otherRead) : -1;
                    }

                    writer.write(String.format(",%s,%d,%d,%s,%d,%d",
                            read.Chromosome, read.PosStart, read.PosEnd, read.Cigar.toString(),
                            read.fragmentInsertSize(), calcFragmentLength));

                    writer.write(String.format(",%s,%d,%s,%s,%s",
                            read.mateChromosome(), read.mateStartPosition(), read.isFirstOfPair(), read.isReadReversed(),
                            read.hasSuppAlignment() ? read.getSuppAlignment() : "NONE"));

                    writer.write(String.format(",%s,%d,%s,%s,%d,%d,%d,%s,%d,%d",
                            geneReadType, transId, transType, validTranscripts,
                            region.getExonRank(transId), region.start(), region.end(), matchType,
                            read.getSoftClipRegionsMatched()[SE_START], read.getSoftClipRegionsMatched()[SE_END]));

                    writer.newLine();
                }
            }
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write read data file: {}", e.toString());
        }
    }

    @VisibleForTesting
    public void processReadRecords(final GeneCollection geneCollection, final List<ReadRecord> readRecords)
    {
        mCurrentGenes = geneCollection;
        mBaseDepth.initialise(geneCollection.regionBounds());

        mValidReadStartRegion[SE_START] = mCurrentGenes.getNonGenicPositions()[SE_START] >= 0
                ? mCurrentGenes.getNonGenicPositions()[SE_START] : mCurrentGenes.regionBounds()[SE_START];

        mValidReadStartRegion[SE_END] = mCurrentGenes.regionBounds()[SE_END];

        mAltSpliceJunctionFinder.setGeneData(mCurrentGenes);
        mRetainedIntronFinder.setGeneData(mCurrentGenes);
        mChimericReads.initialise(mCurrentGenes);

        if(readRecords.size() == 2)
        {
            readRecords.get(0).setFlag(SAMFlag.FIRST_OF_PAIR, true);
        }

        readRecords.forEach(x -> processRead(x));
    }


}
