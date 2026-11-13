package com.hartwig.hmftools.isofox;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.firstInPair;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_PAIR;
import static com.hartwig.hmftools.common.region.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxFunction.ALT_SPLICE_JUNCTIONS;
import static com.hartwig.hmftools.isofox.WriteType.SPLICE_SITE;
import static com.hartwig.hmftools.isofox.common.FragmentMatchType.DISCORDANT;
import static com.hartwig.hmftools.isofox.common.FragmentType.ALT;
import static com.hartwig.hmftools.isofox.common.FragmentType.CHIMERIC;
import static com.hartwig.hmftools.isofox.common.FragmentType.DUPLICATE;
import static com.hartwig.hmftools.isofox.common.FragmentType.FORWARD_STRAND;
import static com.hartwig.hmftools.isofox.common.FragmentType.REVERSE_STRAND;
import static com.hartwig.hmftools.isofox.common.FragmentType.TOTAL;
import static com.hartwig.hmftools.isofox.common.FragmentType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.isofox.common.FragmentType.UNSPLICED;
import static com.hartwig.hmftools.isofox.IsofoxFunction.FUSIONS;
import static com.hartwig.hmftools.isofox.common.Read.findOverlappingRegions;
import static com.hartwig.hmftools.isofox.common.ReadTranscriptUtils.markRegionBases;
import static com.hartwig.hmftools.isofox.common.ReadUtils.consensusDuplicateCount;
import static com.hartwig.hmftools.isofox.common.TransMatchType.SPLICE_JUNCTION;
import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_END;

import static com.hartwig.hmftools.common.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.writeReadData;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.common.BaseDepth;
import com.hartwig.hmftools.isofox.common.Fragment;
import com.hartwig.hmftools.isofox.common.FragmentMatchType;
import com.hartwig.hmftools.isofox.common.FragmentTracker;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.FragmentType;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.Read;
import com.hartwig.hmftools.isofox.common.ReadTranscriptUtils;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.common.TransExonRef;
import com.hartwig.hmftools.isofox.expression.CategoryCountsData;
import com.hartwig.hmftools.isofox.adjusts.GcRatioCounts;
import com.hartwig.hmftools.isofox.expression.ExpressionReadTracker;
import com.hartwig.hmftools.isofox.fusion.ChimericReadTracker;
import com.hartwig.hmftools.isofox.fusion.ChimericUtils;
import com.hartwig.hmftools.isofox.novel.AltSjCohortCache;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder;
import com.hartwig.hmftools.isofox.novel.RetainedIntronFinder;
import com.hartwig.hmftools.isofox.novel.SpliceSiteCounter;
import com.hartwig.hmftools.isofox.results.ResultsWriter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FragmentAllocator
{
    private final IsofoxConfig mConfig;
    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    // state relating to the current gene
    private GeneCollection mCurrentGenes;
    private final FragmentTracker mFragmentReads; // cache of single read until both are available - ie a fragment

    private int mGeneReadCount;
    private int mTotalBamReadCount;
    private int mNextGeneCountLog;

    private final ExpressionReadTracker mExpressionReadTracker;
    private final AltSpliceJunctionFinder mAltSpliceJunctionFinder;
    private final RetainedIntronFinder mRetainedIntronFinder;
    private final ChimericReadTracker mChimericReads;
    private final SpliceSiteCounter mSpliceSiteCounter;
    private final int[] mValidReadStartRegion;
    private final BaseDepth mBaseDepth;

    private final boolean mRunFusions;
    private final boolean mFusionsOnly;
    private final boolean mStatsOnly;

    private final BufferedWriter mReadDataWriter;
    private final BufferedWriter mMultiMapLociWriter;

    private final EnsemblDataCache mGeneTransCache;

    private static final int GENE_LOG_COUNT = 2000000;
    private static final int NON_GENIC_BASE_DEPTH_WIDTH = 250000;

    public FragmentAllocator(
            final IsofoxConfig config, final EnsemblDataCache geneTransCache, final AltSjCohortCache altSjCohortCache,
            final ResultsWriter resultsWriter)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mCurrentGenes = null;
        mFragmentReads = new FragmentTracker();

        mRunFusions = mConfig.Functions.contains(FUSIONS);
        mFusionsOnly = mConfig.runFusionsOnly();
        mStatsOnly = mConfig.runStatisticsOnly();

        mGeneReadCount = 0;
        mTotalBamReadCount = 0;
        mNextGeneCountLog = 0;
        mValidReadStartRegion = new int[SE_PAIR];

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        // reads with supplementary alignment data are only used for chimeric read handling (eg fusions & alt-SJs)
        boolean keepSupplementaries = mRunFusions || mConfig.runFunction(ALT_SPLICE_JUNCTIONS);

        // fusions typically aren't run without expression, but for STAR the existing logic was to drop reads with map qual less than the max
        int minMapQuality = 0;

        // no need to process duplicates since their counts can be extracted from consensus reads for transcript counting
        mBamSlicer = new BamSlicer(minMapQuality, false, keepSupplementaries, false);

        mReadDataWriter = resultsWriter.getReadDataWriter();
        mMultiMapLociWriter = resultsWriter.getMultiMapLociWriter();
        mBaseDepth = new BaseDepth();

        mChimericReads = new ChimericReadTracker(mConfig);
        mChimericReads.setChimericPosDataWriter(resultsWriter.getChimericPositionDataWriter());

        mSpliceSiteCounter = new SpliceSiteCounter(resultsWriter.getSpliceSiteWriter());

        mExpressionReadTracker = new ExpressionReadTracker(mConfig);

        mAltSpliceJunctionFinder = new AltSpliceJunctionFinder(
                mConfig, altSjCohortCache, resultsWriter.getAltSjUnfilteredWriter(), resultsWriter.getAltSjPassingWriter());

        mRetainedIntronFinder = new RetainedIntronFinder(mConfig, resultsWriter.getRetainedIntronWriter());
    }

    public int totalReadCount() { return mTotalBamReadCount; }

    public final GcRatioCounts getGcRatioCounts() { return mExpressionReadTracker.getGcRatioCounts(); }
    public final GcRatioCounts getGeneGcRatioCounts() { return mExpressionReadTracker.getGeneGcRatioCounts(); }
    public BaseDepth getBaseDepth() { return mBaseDepth; }
    public final ChimericReadTracker getChimericReadTracker() { return mChimericReads; }
    public final SpliceSiteCounter getSpliceSiteCounter() { return mSpliceSiteCounter; }

    public void clearCache()
    {
        mFragmentReads.clear();
        mChimericReads.clearData();

        mExpressionReadTracker.setGeneData(null);
        mAltSpliceJunctionFinder.setGeneData(null);
        mRetainedIntronFinder.setGeneData(null);
        mSpliceSiteCounter.clear();

        mCurrentGenes = null;
    }

    public void processBam(final GeneCollection geneCollection, final ChrBaseRegion geneRegion)
    {
        clearCache();

        mCurrentGenes = geneCollection;

        mGeneReadCount = 0;
        mNextGeneCountLog = GENE_LOG_COUNT;

        // and width around the base depth region to pick up junctions outside the gene
        int[] baseDepthRange = new int[SE_PAIR];
        baseDepthRange[SE_START] = max(geneRegion.start(), geneCollection.regionBounds()[SE_START] - NON_GENIC_BASE_DEPTH_WIDTH);
        baseDepthRange[SE_END] = min(geneRegion.end(), geneCollection.regionBounds()[SE_END] + NON_GENIC_BASE_DEPTH_WIDTH);
        mBaseDepth.initialise(baseDepthRange);

        mChimericReads.initialise(mCurrentGenes);

        if(mExpressionReadTracker.enabled())
            mExpressionReadTracker.setGeneData(mCurrentGenes);

        if(mAltSpliceJunctionFinder.enabled())
            mAltSpliceJunctionFinder.setGeneData(mCurrentGenes);

        if(mRetainedIntronFinder.enabled())
            mRetainedIntronFinder.setGeneData(mCurrentGenes);

        mValidReadStartRegion[SE_START] = geneRegion.start();
        mValidReadStartRegion[SE_END] = geneRegion.end();

        mBamSlicer.slice(mSamReader, geneRegion, this::processSamRecord);

        if(mChimericReads.enabled())
            processIncompleteReads();

        ISF_LOGGER.trace("genes({}) bamReadCount({}) depth(bases={} perc={} max={})",
                mCurrentGenes.geneNames(), mGeneReadCount, mBaseDepth.basesWithDepth(),
                String.format("%.3f", mBaseDepth.basesWithDepthPerc()), mBaseDepth.maxDepth());
    }

    private void processSamRecord(final SAMRecord record)
    {
        // to avoid double-processing of reads overlapping 2 (or more) gene collections, only process them if they start in this
        // gene collection or its preceding non-genic region
        if(!positionWithin(record.getStart(), mValidReadStartRegion[SE_START], mValidReadStartRegion[SE_END]))
            return;

        if(mConfig.LogReadIds.contains(record.getReadName()))
        {
            ISF_LOGGER.debug("specific read: {}", readToString(record));
        }

        if(mConfig.Filters.skipRead(record, true))
        {
            ++mChimericReads.getStats().Excluded;
            return;
        }

        trackFragmentCounts(record);

        Read read = new Read(record);

        processRead(read);
    }

    private void trackFragmentCounts(final SAMRecord record)
    {
        ++mTotalBamReadCount;
        ++mGeneReadCount;

        // count each fragment once by only taking the first read, and supplmentaries are ignored
        if(record.getSupplementaryAlignmentFlag())
            return;

        if(!firstInPair(record))
            return;

        mCurrentGenes.addCount(TOTAL, 1);

        if(record.hasAttribute(CONSENSUS_READ_ATTRIBUTE))
        {
            int duplicateCount = consensusDuplicateCount(record);
            mCurrentGenes.addCount(TOTAL, duplicateCount);
            mCurrentGenes.addCount(DUPLICATE, duplicateCount);
        }
    }

    private void processRead(final Read read)
    {
        // for each record find all exons with an overlap
        // skip records if either end isn't in one of the exons for this gene

        if(mGeneReadCount >= mNextGeneCountLog)
        {
            mNextGeneCountLog += GENE_LOG_COUNT;
            ISF_LOGGER.info("chr({}) genes({}) bamRecordCount({})", mCurrentGenes.chromosome(), mCurrentGenes.geneNames(), mGeneReadCount);
        }

        if(reachedGeneReadLimit())
            return;

        mCurrentGenes.setReadGeneCollections(read, mValidReadStartRegion);

        if(read.isSupplementaryAlignment())
        {
            handleSupplementaryRead(read);
        }
        else
        {
            checkFragmentRead(read);
        }
    }

    private boolean checkFragmentRead(final Read read)
    {
        if(!read.isReadPaired())
        {
            processFragmentReads(new Fragment(read));
            return true;
        }

        // check if the 2 reads from a fragment exist and if so handle them a pair, returning true
        Read otherRead = mFragmentReads.checkRead(read);

        if(otherRead != null)
        {
            processFragmentReads(new Fragment(read, otherRead));
            return true;
        }

        return false;
    }

    private void markGeneDataRegions(final Read read)
    {
        List<RegionReadData> overlappingRegions = findOverlappingRegions(mCurrentGenes.getExonRegions(), read);

        if(!overlappingRegions.isEmpty())
        {
            ReadTranscriptUtils.processOverlappingRegions(read, overlappingRegions);
        }
    }

    private void handleSupplementaryRead(final Read read)
    {
        if(read.isMateUnmapped())
            return;

        mBaseDepth.processRead(read.getMappedRegionCoords());

        markGeneDataRegions(read);
        mChimericReads.addSupplementaryRead(read);

        writeChimericReadData(read);
    }

    private void processFragmentReads(final Fragment fragment)
    {
        /* process fragment:
            - fully outside the gene (due to the buffer used, ignore
            - read through a gene ie start or end outside
            - purely intronic
            - chimeric an inversion or translocation
            - supporting 1 or more transcripts
                - fragment fully with an exon - if exon has only 1 transcript then consider unambiguous
                - fragment within 2 exons (including spanning intermediary ones) and/or either exon at the boundary
            - not supporting any transcript - eg alternative splice sites or unspliced reads
        */

        fragment.trimAdapterBases();

        for(Read read : fragment.reads())
        {
            markGeneDataRegions(read);
        }

        int fragmentCount = fragment.fragmentCount();
        int numLoci = fragment.minNumLoci();
        boolean isMultiMapped = numLoci > 1;

        List<Read.AltAlignment> altLoci = fragment.altLoci();

        boolean isChimeric = mChimericReads.isChimeric(fragment, isMultiMapped);

        if(mStatsOnly)
        {
            if(isChimeric)
                mCurrentGenes.addCount(CHIMERIC, 1);
            else if(fragment.isFullyIntronic())
                mCurrentGenes.addCount(UNSPLICED, fragmentCount);
            else
                mCurrentGenes.addCount(TRANS_SUPPORTING, fragmentCount);

            return;
        }

        List<BaseRegion> commonMappings = fragment.mergedMappings();

        mBaseDepth.processRead(commonMappings);

        // if either read is chimeric (including one outside the genic region) then handle them both as such
        // some of these may be re-processed as alternative SJ candidates if they are within a single gene
        if(isChimeric)
        {
            if(!isMultiMapped)
            {
                if(mChimericReads.enabled())
                    mChimericReads.addChimericFragment(fragment);
                else
                    mCurrentGenes.addCount(CHIMERIC, 1);
            }

            if(mReadDataWriter != null && mConfig.writeType(WriteType.READ))
            {
                List<GeneReadData> overlapGenes = mCurrentGenes.findGenesCoveringRange(
                        fragment.minAlignmentStart(),
                        fragment.maxAlignmentEnd(), true);

                for(Read read : fragment.reads())
                {
                    writeReadData(mReadDataWriter, overlapGenes, read, CHIMERIC, 0);
                }
            }

            return;
        }

        if(mRunFusions)
        {
            // reads with sufficient soft-clipping and not mapped to an adjacent region are candidates for fusion re-alignment
            if(fragment.reads().stream().anyMatch(ChimericUtils::isRealignedFragmentCandidate))
            {
                mChimericReads.addRealignmentCandidates(fragment);
            }

            if(mFusionsOnly)
                return;
        }

        if(numLoci > 1 && altLoci != null && mMultiMapLociWriter != null)
            recordMultiMapLoci(fragment, altLoci);

        int readPosMin = fragment.minAlignmentStart();
        int readPosMax = fragment.maxAlignmentEnd();

        List<GeneReadData> overlapGenes = mCurrentGenes.findGenesCoveringRange(readPosMin, readPosMax, true);

        if(fragment.isFullyIntronic())
        {
            // fully intronic read in every transcript and gene
            processIntronicReads(overlapGenes, fragment, fragmentCount, isMultiMapped);
            return;
        }

        // first find valid transcripts in both reads
        List<Integer> validTranscripts = Lists.newArrayList();

        List<RegionReadData> validRegions = fragment.uniqueValidRegions();

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
        if(mConfig.writeType(SPLICE_SITE))
        {
            mSpliceSiteCounter.registerSpliceSiteSupport(fragment, mCurrentGenes.getExonRegions());
        }

        for(int transId : fragment.validTypeTranscripts())
        {
            int calcFragmentLength = calcFragmentLength(transId, fragment);

            if(calcFragmentLength > 0 && calcFragmentLength <= mConfig.MaxFragmentLength)
                validTranscripts.add(transId);
        }

        Set<Integer> invalidTranscripts = fragment.invalidTranscripts(validTranscripts);

        FragmentType fragmentType = UNSPLICED;

        // now mark all other transcripts which aren't valid either due to the read pair
        if(validTranscripts.isEmpty())
        {
            // no valid transcripts but record against the gene further information about these reads
            boolean checkRetainedIntrons = false;

            if(fragment.containsSplit())
            {
                fragmentType = ALT;

                if(mAltSpliceJunctionFinder.enabled())
                {
                    mAltSpliceJunctionFinder.evaluateFragmentReads(
                            overlapGenes, fragment.reads(), invalidTranscripts.stream().collect(Collectors.toList()));
                }

                checkRetainedIntrons = true;
            }
            else
            {
                // look for alternative splicing from long reads involving more than one region and not spanning into an intron
                for(int transId : invalidTranscripts)
                {
                    if(fragment.spansMultipleRegions(transId))
                    {
                        fragmentType = ALT;
                        break;
                    }
                }

                checkRetainedIntrons = true;
            }

            if(checkRetainedIntrons && mRetainedIntronFinder.enabled())
                mRetainedIntronFinder.evaluateFragmentReads(fragment);

            if(fragmentType == UNSPLICED)
            {
                mExpressionReadTracker.processUnsplicedGenes(overlapGenes, validTranscripts, commonMappings, fragmentCount, isMultiMapped);
            }
        }
        else
        {
            // record valid read info against each region now that it is known
            fragmentType = TRANS_SUPPORTING;

            // first mark any invalid trans as 'other' meaning it doesn't require any further classification since a valid trans exists
            fragment.setOtherTranscripts(validTranscripts);

            if(mConfig.RunValidations)
            {
                for(BaseRegion readRegion : commonMappings)
                {
                    if(commonMappings.stream().filter(x -> x.start() == readRegion.start() && x.end() == readRegion.end()).count() > 1)
                    {
                        ISF_LOGGER.error("repeated read region({})", readRegion);
                    }
                }
            }

            validRegions.forEach(x -> markRegionBases(commonMappings, x));

            // now set counts for each valid transcript
            boolean isUniqueTrans = validTranscripts.size() == 1;
            Boolean supportedGeneIsForward = null;

            FragmentMatchType comboTransMatchType = FragmentMatchType.SHORT;

            for(int transId : validTranscripts)
            {
                int regionCount = (int)validRegions.stream().filter(x -> x.hasTransId(transId)).count();

                FragmentMatchType transMatchType;

                if(supportedGeneIsForward == null)
                {
                    for(Read read : fragment.reads())
                    {
                        supportedGeneIsForward = findGeneStrand(read, validTranscripts);

                        if (supportedGeneIsForward != null)
                            break;
                    }
                }

                if(fragment.hasTranscriptClassification(transId, SPLICE_JUNCTION))
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
                if(!fragment.containsSplit() && fragment.readsInDifferentExons(transId))
                {
                    mCurrentGenes.addTranscriptReadMatch(transId, DISCORDANT);
                }

                // keep track of which regions have been allocated from this fragment as a whole, so not counting each read separately
                mExpressionReadTracker.processValidTranscript(transId, fragment.reads(), isUniqueTrans);
            }

            mExpressionReadTracker.processUnsplicedGenes(
                    comboTransMatchType, overlapGenes, validTranscripts, commonMappings, fragmentCount, isMultiMapped);

            if(supportedGeneIsForward != null)
            {
                // track fragment strandedness
                Orientation fragmentOrientation = fragment.orientation();

                if(fragmentOrientation != null)
                {
                    if(fragmentOrientation.isForward() == supportedGeneIsForward)
                        mCurrentGenes.addCount(FORWARD_STRAND, fragmentCount);
                    else
                        mCurrentGenes.addCount(REVERSE_STRAND, fragmentCount);
                }
            }
        }

        mCurrentGenes.addCount(fragmentType, fragmentCount);

        if(mReadDataWriter != null && mConfig.writeType(WriteType.READ))
        {
            for(Read read : fragment.reads())
            {
                writeReadData(mReadDataWriter, overlapGenes, read, fragmentType, validTranscripts.size());
            }
        }
    }

    private Boolean findGeneStrand(final Read read, final List<Integer> transcripts)
    {
        for(int transId : transcripts)
        {
            for(RegionReadData regionReadData : read.getMappedRegions().keySet())
            {
                TransExonRef transExonRef = regionReadData.getTransExonRefs().stream().filter(x -> x.TransId == transId).findFirst().orElse(null);
                if(transExonRef != null)
                {
                    GeneReadData geneData = mCurrentGenes.genes().stream()
                            .filter(x -> x.Gene.GeneId.equals(transExonRef.GeneId)).findFirst().orElse(null);

                    if(geneData != null)
                        return geneData.Gene.Strand == ORIENT_FWD;
                }
            }
        }

        return null;
    }

    private int calcFragmentLength(int transId, final Fragment fragment)
    {
        TranscriptData transData = mCurrentGenes.getTranscripts().stream().filter(x -> x.TransId == transId).findFirst().orElse(null);
        if(transData == null)
            return -1;

        return ReadTranscriptUtils.calcFragmentLength(transData, fragment.minAlignmentStart(), fragment.maxAlignmentEnd());
    }

    private boolean reachedGeneReadLimit()
    {
        if(mConfig.GeneReadLimit == 0 || mGeneReadCount < mConfig.GeneReadLimit)
            return false;

        mBamSlicer.haltProcessing();
        ISF_LOGGER.warn("chr({}) genes({}) readCount({}) exceeds max read count",
                mCurrentGenes.chromosome(), mCurrentGenes.geneNames(), mGeneReadCount);
        return true;
    }

    public List<CategoryCountsData> getTransComboData() { return mExpressionReadTracker.getTransComboData(); }

    private boolean altOverlapsExon(final GeneData gene, final ChrBaseRegion altRegion)
    {
        List<TranscriptData> transcripts = mGeneTransCache.getTranscripts(gene.GeneId);

        if(transcripts == null)
            return false;

        for(TranscriptData transData : transcripts)
        {
            for(ExonData exon : transData.exons())
            {
                if(positionsOverlap(exon.Start, exon.End, altRegion.start(), altRegion.end()))
                    return true;
            }
        }

        return false;
    }

    private void recordMultiMapLoci(final Fragment fragment, final List<Read.AltAlignment> altLoci)
    {
        if(mMultiMapLociWriter == null)
            return;

        // record a multi-mapped fragment's primary alignment plus each XA alternate locus (opt-in WriteType.MULTI_MAP_LOCI);
        // InGeneCollection flags whether the locus falls within the gene collection currently being processed
        int fragStart = fragment.minAlignmentStart();
        int fragEnd = fragment.maxAlignmentEnd();
        boolean primarySpliced = fragment.containsSplit();

        writeMultiMapLocus(
                mMultiMapLociWriter, mCurrentGenes.id(), fragment.id(), "PRIMARY", fragment.chromosome(), fragStart, fragEnd,
                primarySpliced, mCurrentGenes.geneNames(), true);

        int[] bounds = mCurrentGenes.regionBounds();

        for(Read.AltAlignment locus : altLoci)
        {
            boolean inGeneCollection = locus.Region.Chromosome.equals(mCurrentGenes.chromosome())
                    && positionsOverlap(locus.Region.start(), locus.Region.end(), bounds[SE_START], bounds[SE_END]);

            writeMultiMapLocus(
                    mMultiMapLociWriter, mCurrentGenes.id(), fragment.id(), "XA", locus.Region.Chromosome,
                    locus.Region.start(), locus.Region.end(), locus.Spliced, altExonicGeneNames(locus), inGeneCollection);
        }
    }

    private String altExonicGeneNames(final Read.AltAlignment locus)
    {
        if(mGeneTransCache == null)
            return "";

        return mGeneTransCache.findGeneByRange(locus.Region.Chromosome, locus.Region.start(), locus.Region.end()).stream()
                .filter(gene -> altOverlapsExon(gene, locus.Region))
                .map(gene -> gene.GeneName)
                .collect(Collectors.joining(ITEM_DELIM));
    }

    private synchronized static void writeMultiMapLocus(
            final BufferedWriter writer, int geneCollectionId, final String readId, final String recordType,
            final String chromosome, int posStart, int posEnd, boolean spliced, final String genes, boolean inGeneCollection)
    {
        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(String.valueOf(geneCollectionId));
            sj.add(readId);
            sj.add(recordType);
            sj.add(chromosome);
            sj.add(String.valueOf(posStart));
            sj.add(String.valueOf(posEnd));
            sj.add(String.valueOf(spliced));
            sj.add(genes);
            sj.add(String.valueOf(inGeneCollection));
            writer.write(sj.toString());
            writer.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write multi-map loci data: {}", e.toString());
        }
    }

    private void processIntronicReads(
            final List<GeneReadData> genes, final Fragment fragment, int fragmentCount, boolean multiMapped)
    {
        if(fragment.containsSplit())
        {
            mCurrentGenes.addCount(ALT, 1); // does not count duplicates since not expression related

            if(mAltSpliceJunctionFinder.enabled())
                mAltSpliceJunctionFinder.evaluateFragmentReads(genes, fragment.reads(), Lists.newArrayList());

            return;
        }

        mExpressionReadTracker.processIntronicReads(genes, fragment, fragmentCount, multiMapped);
        mCurrentGenes.addCount(UNSPLICED, fragmentCount);

        if(mReadDataWriter != null && mConfig.writeType(WriteType.READ))
        {
            for(Read read  : fragment.reads())
            {
                writeReadData(mReadDataWriter, genes, read, UNSPLICED, 0);
            }
        }
    }

    private void processIncompleteReads()
    {
        // now slicing has completed, process unpaired primaries (supps have been handled already)
        for(Object readObject : mFragmentReads.readMap().values())
        {
            Read read = (Read)readObject;
            markGeneDataRegions(read);
            mBaseDepth.processRead(read.getMappedRegionCoords());
            writeChimericReadData(read);
        }

        mChimericReads.postProcessChimericReads(mBaseDepth, mFragmentReads);
        processChimericNovelJunctions();
    }

    private void processChimericNovelJunctions()
    {
        // examine chimeric reads to see if they can instead be handled as novel alternate splicing
        if(!mAltSpliceJunctionFinder.enabled() || mChimericReads.getLocalChimericReads().isEmpty())
            return;

        List<Integer> invalidTrans = Lists.newArrayList();

        for(List<Read> reads : mChimericReads.getLocalChimericReads())
        {
            Read read1 = null;
            Read read2 = null;

            if(reads.size() == 2)
            {
                read1 = reads.get(0);
                read2 = reads.get(1);
            }
            else if(reads.size() == 3)
            {
                for(Read read : reads)
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

            int readPosMin = min(read1.alignmentStart(), read2.alignmentStart());
            int readPosMax = max(read1.alignmentEnd(), read2.alignmentEnd());

            List<GeneReadData> overlapGenes = mCurrentGenes.findGenesCoveringRange(readPosMin, readPosMax, false);
            mAltSpliceJunctionFinder.evaluateFragmentReads(overlapGenes, List.of(read1, read2), invalidTrans);
        }
    }

    public void annotateNovelLocations()
    {
        recordNovelLocationReadDepth();

        if(mAltSpliceJunctionFinder.enabled())
        {
            mAltSpliceJunctionFinder.prioritiseGenes();
            mAltSpliceJunctionFinder.finalise();
            mAltSpliceJunctionFinder.writeAltSpliceJunctions();
        }

        if(mRetainedIntronFinder.enabled())
            mRetainedIntronFinder.writeRetainedIntrons();
    }

    private void recordNovelLocationReadDepth()
    {
        if(mAltSpliceJunctionFinder.getAltSpliceJunctions().isEmpty() && mRetainedIntronFinder.getRetainedIntrons().isEmpty())
            return;

        if(mAltSpliceJunctionFinder.enabled())
            mAltSpliceJunctionFinder.setPositionDepth(mBaseDepth);

        if(mRetainedIntronFinder.enabled())
            mRetainedIntronFinder.setPositionDepth(mBaseDepth);
    }

    public void registerKnownFusionPairs(final EnsemblDataCache geneTransCache)
    {
        mChimericReads.registerKnownFusionPairs(geneTransCache);
    }

    private void writeChimericReadData(final Read read)
    {
        if(mReadDataWriter == null)
            return;

        List<GeneReadData> overlapGenes = mCurrentGenes.findGenesCoveringRange(
                read.alignmentStart(), read.alignmentEnd(), true);

        writeReadData(mReadDataWriter, overlapGenes, read, CHIMERIC, 0);
    }

    @VisibleForTesting
    public void processReadRecords(final GeneCollection geneCollection, final List<Read> reads)
    {
        mCurrentGenes = geneCollection;
        mBaseDepth.initialise(geneCollection.regionBounds());

        mValidReadStartRegion[SE_START] = mCurrentGenes.getNonGenicPositions()[SE_START] >= 0
                ? mCurrentGenes.getNonGenicPositions()[SE_START] : mCurrentGenes.regionBounds()[SE_START];

        mValidReadStartRegion[SE_END] = mCurrentGenes.regionBounds()[SE_END];

        mExpressionReadTracker.setGeneData(mCurrentGenes);
        mAltSpliceJunctionFinder.setGeneData(mCurrentGenes);
        mRetainedIntronFinder.setGeneData(mCurrentGenes);
        mChimericReads.initialise(mCurrentGenes);

        reads.forEach(x -> processRead(x));
    }

    @VisibleForTesting
    public void postSliceProcessReads()
    {
        processIncompleteReads();
    }

    @VisibleForTesting
    public final FragmentTracker getFragmentTracker() { return mFragmentReads; }

}
