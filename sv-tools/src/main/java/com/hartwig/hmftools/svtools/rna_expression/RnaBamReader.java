package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.CHIMERIC;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.DUPLICATE;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.READ_THROUGH;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.TOTAL;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.UNSPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_LONG;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_SHORT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.TC_SPLICE;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.getUniqueValidRegion;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.hasSkippedExons;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.markRegionBases;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.validRegionMatchType;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.validTranscriptType;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.INTRONIC;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.deriveCommonRegions;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.ALT;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.OTHER_TRANS;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.SPLICE_JUNCTION;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.SAMSlicer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class RnaBamReader
{
    private final RnaExpConfig mConfig;
    private final SamReader mSamReader;

    private FragmentSizeCalcs mFragmentSizeCalcs;

    // state relating to the current gene
    private GeneReadData mCurrentGene;
    private final List<String> mDiscardedReads;

    private int mGeneReadCount;
    private int mTotalBamReadCount;

    private static final int DEFAULT_MIN_MAPPING_QUALITY = 1;

    private final Map<String,ReadRecord> mFragmentReads;
    private final Map<Long,List<long[]>> mDuplicateCache;
    private final List<String> mDuplicateReadIds;

    private BufferedWriter mWriter;

    private static final Logger LOGGER = LogManager.getLogger(RnaBamReader.class);

    public RnaBamReader(final RnaExpConfig config)
    {
        mConfig = config;

        mCurrentGene = null;
        mFragmentReads = Maps.newHashMap();
        mDiscardedReads = Lists.newArrayList();

        mGeneReadCount = 0;
        mTotalBamReadCount = 0;

        if(!mConfig.BamFile.isEmpty() && Files.exists(Paths.get(mConfig.BamFile)))
        {
            mSamReader = SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile));
        }
        else
        {
            mSamReader = null;
        }

        mDuplicateCache = Maps.newHashMap();
        mDuplicateReadIds = Lists.newArrayList();

        mWriter = null;
    }

    public boolean validReader() { return mSamReader != null; }

    public void setFragmentSizeCalcs(final FragmentSizeCalcs calcs) { mFragmentSizeCalcs = calcs; }

    public void close()
    {
        LOGGER.info("read {} total BAM records", mTotalBamReadCount);

        closeBufferedWriter(mWriter);
    }

    public void readBamCounts(final GeneReadData geneReadData, final GenomeRegion genomeRegion)
    {
        mFragmentReads.clear();
        mDiscardedReads.clear();

        clearDuplicates();

        mCurrentGene = geneReadData;
        mGeneReadCount = 0;

        SAMSlicer samSlicer = new SAMSlicer(DEFAULT_MIN_MAPPING_QUALITY, Lists.newArrayList(genomeRegion));
        samSlicer.setDropDuplicates(false);

        samSlicer.slice(mSamReader, this::processSamRecord);

        if(!mFragmentReads.isEmpty())
        {
            if(LOGGER.isDebugEnabled())
            {
                int discardedReadMatches = (int) mFragmentReads.keySet().stream().filter(x -> mDiscardedReads.contains(x)).count();

                LOGGER.debug("gene({}) has {} unmatched reads, discarded({})",
                        mCurrentGene.GeneData.GeneName, mFragmentReads.size(), discardedReadMatches);
            }

            // only paired reads can be processed - otherwise implies a link to another gene or non-genic region

            // mFragmentReads.values().forEach(x -> processingUnmatchedRead(x));
            mFragmentReads.clear();
        }

        LOGGER.debug("gene({}) bamReadCount({})", mCurrentGene.GeneData.GeneName, mGeneReadCount);

        if(mConfig.GeneStatsOnly)
        {
            mCurrentGene.addCount(TOTAL, mGeneReadCount / 2);
        }
    }

    public void readBamCounts(final GenomeRegion genomeRegion, final Consumer<SAMRecord> consumer)
    {
        clearDuplicates();
        SAMSlicer samSlicer = new SAMSlicer(DEFAULT_MIN_MAPPING_QUALITY, Lists.newArrayList(genomeRegion));
        samSlicer.slice(mSamReader, consumer);
    }

    private void processSamRecord(@NotNull final SAMRecord record)
    {
        if(checkDuplicates(record))
        {
            if(record.getFirstOfPairFlag())
                mCurrentGene.addCount(DUPLICATE, 1);

            if(!mConfig.KeepDuplicates)
                return;
        }

        ++mTotalBamReadCount;
        ++mGeneReadCount;

        // recordFragmentLength(record); // now an stand-alone routine

        if(mConfig.GeneStatsOnly)
            return;

        processRead(ReadRecord.from(record));
    }

    public void processRead(ReadRecord read)
    {
        // for each record find all exons with an overlap
        // skip records if either end isn't in one of the exons for this gene

        if(mGeneReadCount > 0 && (mGeneReadCount % 100000) == 0)
        {
            LOGGER.debug("gene({}) bamRecordCount({})", mCurrentGene.GeneData.GeneName, mGeneReadCount);
        }

        if(mConfig.ReadCountLimit > 0 && mGeneReadCount >= mConfig.ReadCountLimit)
        {
            if(mGeneReadCount == mConfig.ReadCountLimit)
            {
                LOGGER.warn("gene({}) readCount({}) exceeds max read count", mCurrentGene.GeneData.GeneName, mGeneReadCount);
            }

            return;
        }

        if(read.translocation())
        {
            mCurrentGene.addCount(TOTAL, 1);
            mCurrentGene.addCount(CHIMERIC, 1);
            return;
        }

        boolean exonOverlap = mCurrentGene.getExonRegions().stream()
                .anyMatch(x -> !(read.PosEnd < x.start() || read.PosStart > x.Region.end()));

        if(!exonOverlap)
        {
            boolean outsideGene = read.PosStart > mCurrentGene.GeneData.GeneEnd || read.PosEnd < mCurrentGene.GeneData.GeneStart;

            if(outsideGene)
            {
                checkFragmentRead(read);
            }
            else
            {
                checkIntronicRegions(read);
            }

            return;
        }

        // the read is fully within the exon
        List<RegionReadData> overlappingRegions = mCurrentGene.getExonRegions().stream()
                .filter(x -> read.overlapsMappedReads(x.Region.start(), x.Region.end()))
                .collect(Collectors.toList());

        if(!overlappingRegions.isEmpty())
        {
            // look at all matched reads within the context of a transcript
            read.processOverlappingRegions(overlappingRegions);
        }

        checkFragmentRead(read);
    }

    private void processFragmentReads(@NotNull final ReadRecord read1, @NotNull final ReadRecord read2)
    {
        /* use of fragment read pair:
            - supporting a transcript:
                - both reads fully with an exon - if exon has only 1 transcript then consider unambiguous
                - both reads within 2 exons (including spanning intermediary ones) and/or either exon at the boundary
            - not supporting a transcript
                - both reads touch the same exon if there is a gap in the reads
                - one read in an intron -> UNSPLICED
        */

        boolean r1OutsideGene = read1.PosStart > mCurrentGene.GeneData.GeneEnd || read1.PosEnd < mCurrentGene.GeneData.GeneStart;
        boolean r2OutsideGene = read2.PosStart > mCurrentGene.GeneData.GeneEnd || read2.PosEnd < mCurrentGene.GeneData.GeneStart;

        if(r1OutsideGene && r2OutsideGene)
            return;

        mCurrentGene.addCount(TOTAL, 1);

        if(read1.localInversion() || read2.localInversion())
        {
            mCurrentGene.addCount(CHIMERIC, 1);
            return;
        }

        if(read1.getMappedRegions().isEmpty() && read2.getMappedRegions().isEmpty())
            return;

        final Map<String,TransMatchType> firstReadTransTypes = read1.getTranscriptClassifications();

        final Map<String,TransMatchType> secondReadTransTypes = read2.getTranscriptClassifications();

        // first find valid transcripts in both reads
        final List<String> validTranscripts = Lists.newArrayList();
        final List<String> invalidTranscripts = Lists.newArrayList();

        final List<RegionReadData> validRegions = getUniqueValidRegion(read1, read2);

        for(Map.Entry<String,TransMatchType> entry : firstReadTransTypes.entrySet())
        {
            final String trans = entry.getKey();

            if(validTranscriptType(entry.getValue()))
            {
                if(secondReadTransTypes.containsKey(trans) && validTranscriptType(secondReadTransTypes.get(trans)))
                {
                    if(!hasSkippedExons(validRegions, trans, mConfig.LongFragmentLimit))
                    {
                        validTranscripts.add(trans);
                        continue;
                    }
                }
            }

            if(!invalidTranscripts.contains(trans))
                invalidTranscripts.add(trans);
        }

        boolean isLongFragment = abs(read1.fragmentInsertSize()) > mConfig.LongFragmentLimit;
        GeneMatchType geneReadType = UNSPLICED;

        // now mark all other transcripts which aren't valid either due to the read pair
        if(validTranscripts.isEmpty())
        {
            // no valid transcripts but record against the gene further information about these reads
            if(r1OutsideGene || r2OutsideGene)
            {
                geneReadType = READ_THROUGH;
            }
            else if(read1.containsSplit() || read2.containsSplit())
            {
                geneReadType = GeneMatchType.ALT;
            }
            else if(isLongFragment)
            {
                // look for alternative splicing from long reads involving more than one region and not spanning into an intron
                for(String trans : invalidTranscripts)
                {
                    List<RegionReadData> regions = read1.getMappedRegions().entrySet().stream()
                            .filter(x -> x.getKey().hasTransId(trans))
                            .filter(x -> x.getValue() != EXON_INTRON)
                            .map(x -> x.getKey()).collect(Collectors.toList());;

                    final List<RegionReadData> regions2 = read2.getMappedRegions().entrySet().stream()
                            .filter(x -> x.getKey().hasTransId(trans))
                            .filter(x -> x.getValue() != EXON_INTRON)
                            .map(x -> x.getKey()).collect(Collectors.toList());

                    for(RegionReadData region : regions2)
                    {
                        if (!regions.contains(region))
                            regions.add(region);
                    }

                    if(regions.size() > 1)
                    {
                        geneReadType = GeneMatchType.ALT;
                        break;
                    }
                }
            }
        }
        else
        {
            // record valid read info against each region now that it is known
            geneReadType = TRANS_SUPPORTING;

            // first mark any invalid trans as 'other' meaning it doesn't require any further classification since a valid trans exists
            firstReadTransTypes.entrySet().stream()
                    .filter(x -> validTranscriptType(x.getValue()))
                    .filter(x -> !validTranscripts.contains(x.getKey()))
                    .forEach(x -> x.setValue(OTHER_TRANS));

            secondReadTransTypes.entrySet().stream()
                    .filter(x -> validTranscriptType(x.getValue()))
                    .filter(x -> !validTranscripts.contains(x.getKey()))
                    .forEach(x -> x.setValue(OTHER_TRANS));

            // now record the bases covered by the read in these matched regions
            final List<long[]> commonMappings = deriveCommonRegions(read1.getMappedRegionCoords(), read2.getMappedRegionCoords());

            validRegions.forEach(x -> markRegionBases(commonMappings, x));

            // now set counts for each valid transcript
            boolean isUniqueTrans = validTranscripts.size() == 1;

            for (final String trans : validTranscripts)
            {
                int regionCount = (int)validRegions.stream().filter(x -> x.hasTransId(trans)).count();

                int transMatchType;

                if(read1.getTranscriptClassification(trans) == SPLICE_JUNCTION || read2.getTranscriptClassification(trans) == SPLICE_JUNCTION)
                    transMatchType = TC_SPLICE;
                else if(regionCount > 1 && isLongFragment)
                    transMatchType = TC_LONG;
                else
                    transMatchType = TC_SHORT;

                mCurrentGene.addTranscriptReadMatch(trans, isUniqueTrans, transMatchType);

                // keep track of which regions have been allocated from this fragment as a whole, so not counting each read separately
                final List<RegionReadData> processedRegions = Lists.newArrayList();

                processValidTranscript(trans, read1, processedRegions, isUniqueTrans);
                processValidTranscript(trans, read2, processedRegions, isUniqueTrans);
            }
        }

        mCurrentGene.addCount(geneReadType, 1);

        if(mConfig.WriteReadData)
        {
            writeReadData(0, read1, geneReadType);
            writeReadData(1, read2, geneReadType);
        }
    }

    private void processValidTranscript(
            final String trans, final ReadRecord read, final List<RegionReadData> processedRegions, boolean isUniqueTrans)
    {
        List<RegionReadData> regions = read.getMappedRegions().entrySet().stream()
                .filter(x -> x.getKey().hasTransId(trans))
                .filter(x -> validRegionMatchType(x.getValue()))
                .map(x -> x.getKey()).collect(Collectors.toList());

        for(RegionReadData region : regions)
        {
            if (!processedRegions.contains(region))
            {
                // register a read against this valid transcript region
                region.addTranscriptReadMatch(trans, isUniqueTrans);
            }
        }

        // any adjacent reads can record a splice junction count
        if(regions.size() > 1 && read.getTranscriptClassification(trans) == SPLICE_JUNCTION)
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
                        region1.addTranscriptJunctionMatch(trans, SE_END, isUniqueTrans);
                        region2.addTranscriptJunctionMatch(trans, SE_START, isUniqueTrans);
                    }
                    else if(region1.getPreRegions().contains(region2))
                    {
                        region1.addTranscriptJunctionMatch(trans, SE_START, isUniqueTrans);
                        region2.addTranscriptJunctionMatch(trans, SE_END, isUniqueTrans);
                    }
                }
            }
        }

        regions.forEach(x -> processedRegions.add(x));
    }

    public static boolean overlaps(final GenomeRegion region, final ReadRecord record)
    {
        // overlapping but neither wholy contained within
        if(region.start() >= record.PosStart && region.start() <= record.PosEnd && region.end() > record.PosEnd)
            return true; // region starts at or within and ends after

        if(region.end() >= record.PosStart && region.end() <= record.PosEnd && region.start() < record.PosStart)
            return true; // region starts before and ends at the record end or before

        return false;
    }

    private void checkIntronicRegions(final ReadRecord read)
    {
        if(read.Cigar == null)
            return;

        if(read.Cigar.containsOperator(CigarOperator.N) || !read.Cigar.containsOperator(CigarOperator.M))
            return;

        /* unused for now
        RegionReadData intronReadData = mCurrentGene.getIntronRegions().stream()
                .filter(x -> read.PosStart >= x.Region.start() && read.PosEnd <= x.Region.end())
                .findFirst().orElse(null);

        if(intronReadData != null)
        {
            if (mConfig.AllTranscripts && intronReadData.getRefRegions().size() == 1)
            {
                // only record intronic reads if they are unique to a transcript
                intronReadData.addMatchedRead(INTRONIC);
            }
        }
        */

        // process the fragment if both reads are now available, and implies one of the reads covers an exon
        if(mFragmentReads.containsKey(read.Id))
        {
            checkFragmentRead(read);
            return;
        }

        // cache this read if it's pair is expected to reach an exon with its pair
        // (for testing assume that the first read encountered is the lower of the 2)
        long otherReadStartPos = read.samRecord() != null ? read.samRecord().getMateAlignmentStart() : read.PosEnd + read.fragmentInsertSize();
        long otherReadEndPos = otherReadStartPos + read.Length; // assume similar length

        // measure distance to nearest exon region and cache if within range of being a fragment read pair
        boolean otherReadExonic = mCurrentGene.getExonRegions().stream()
                .anyMatch(x -> (otherReadStartPos >= x.start() && otherReadStartPos <= x.end())
                        || (otherReadEndPos >= x.start() && otherReadEndPos <= x.end()));

        if(otherReadExonic)
        {
            // cache the read until the exonic-read is processed
            checkFragmentRead(read);
            return;
        }

        if(read.PosStart < otherReadStartPos)
        {
            mCurrentGene.addCount(UNSPLICED, 1);
            mCurrentGene.addCount(TOTAL, 1);
        }

        if(LOGGER.isDebugEnabled())
        {
            if(mDiscardedReads.contains(read.Id))
            {
                // both reads intronic so ignore
                mDiscardedReads.remove(read.Id);
                mFragmentReads.remove(read.Id);
            }
            else
            {
                mDiscardedReads.add(read.Id);
            }
        }
    }

    private boolean checkFragmentRead(ReadRecord read)
    {
        // check if the 2 reads from a fragment exist and if so handle them a pair, returning true
        if(read.samRecord() != null)
        {
            if(!read.samRecord().getMateReferenceName().equals(read.Chromosome)
            || read.samRecord().getMateReferenceIndex() == null)
            {
                return false;
            }
        }

        ReadRecord otherRead = mFragmentReads.get(read.Id);

        if(otherRead != null)
        {
            mFragmentReads.remove(read.Id);
            processFragmentReads(read, otherRead);
            return true;
        }

        mFragmentReads.put(read.Id, read);
        return false;
    }

    private static final int DUP_DATA_SECOND_START = 0;
    private static final int DUP_DATA_READ_LEN = 1;
    private static final int DUP_DATA_INSERT_SIZE = 2;

    public boolean checkDuplicates(final SAMRecord record)
    {
        if(record.getDuplicateReadFlag())
            return true;

        if(!mConfig.MarkDuplicates)
            return false;

        if(mDuplicateReadIds.contains(record.getReadName()))
        {
            mDuplicateReadIds.remove(record.getReadName());
            return true;
        }

        if(!record.getReferenceName().equals(record.getMateReferenceName()) || record.getReadNegativeStrandFlag() == record.getMateNegativeStrandFlag())
            return false;

        long firstStartPos = record.getFirstOfPairFlag() ? record.getStart() : record.getMateAlignmentStart();
        long secondStartPos = record.getFirstOfPairFlag() ? record.getMateAlignmentStart() : record.getStart();
        int readLength = record.getReadLength();
        int insertSize = record.getInferredInsertSize();

        List<long[]> dupDataList = mDuplicateCache.get(firstStartPos);

        if(dupDataList == null)
        {
            dupDataList = Lists.newArrayList();
            mDuplicateCache.put(firstStartPos, dupDataList);
        }
        else
        {
            // search for a match
            if(dupDataList.stream().anyMatch(x -> x[DUP_DATA_SECOND_START] == secondStartPos
                    && x[DUP_DATA_READ_LEN] == readLength && insertSize == x[DUP_DATA_INSERT_SIZE]))
            {
                LOGGER.trace("duplicate fragment: id({}) chr({}) pos({}->{}) otherReadStart({}) insertSize({})",
                        record.getReadName(), record.getReferenceName(), firstStartPos, record.getEnd(), secondStartPos, insertSize);

                // cache so the second read can be identified immediately
                mDuplicateReadIds.add(record.getReadName());
                return true;
            }
        }

        long[] dupData = {secondStartPos, readLength, insertSize};
        dupDataList.add(dupData);

        return false;
    }

    private void clearDuplicates()
    {
        mDuplicateCache.clear();
        mDuplicateReadIds.clear();
    }

    private void recordFragmentLength(final SAMRecord record)
    {
        if(mFragmentSizeCalcs != null && mFragmentSizeCalcs.enabled())
            mFragmentSizeCalcs.recordFragmentLength(record, mCurrentGene);
    }

    private void writeReadData(int readIndex, final ReadRecord read, GeneMatchType geneReadType)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mWriter == null)
            {
                final String outputFileName = mConfig.OutputDir + "RNA_READ_DATA.csv";

                mWriter = createBufferedWriter(outputFileName, false);
                mWriter.write("GeneId,GeneName,ReadIndex,ReadId,Chromosome,PosStart,PosEnd,Cigar");
                mWriter.write(",GeneClass,TransId,TransClass,ExonRank,ExonStart,ExonEnd,RegionClass");
                mWriter.newLine();
            }

            for(Map.Entry<String,TransMatchType> entry : read.getTranscriptClassifications().entrySet())
            {
                final String trans = entry.getKey();
                TransMatchType transType = entry.getValue();

                for(Map.Entry<RegionReadData,RegionMatchType> rEntry : read.getMappedRegions().entrySet())
                {
                    RegionReadData region = rEntry.getKey();
                    RegionMatchType matchType = rEntry.getValue();

                    if(!region.hasTransId(trans))
                        continue;

                    mWriter.write(String.format("%s,%s,%d,%s,%s,%d,%d,%s",
                            mCurrentGene.GeneData.GeneId, mCurrentGene.GeneData.GeneName, readIndex, read.Id,
                            read.Chromosome, read.PosStart, read.PosEnd, read.Cigar.toString()));

                    mWriter.write(String.format(",%s,%s,%s,%d,%d,%d,%s",
                            geneReadType, trans, transType, region.getExonRank(trans), region.start(), region.end(), matchType));

                    mWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            LOGGER.error("failed to write read data file: {}", e.toString());
        }
    }

    @VisibleForTesting
    public void processReadRecords(final GeneReadData geneReadData, final List<ReadRecord> readRecords)
    {
        mCurrentGene = geneReadData;

        readRecords.forEach(x -> processRead(x));
    }


}
