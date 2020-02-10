package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.abs;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_ALT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_CHIMERIC;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_DUPLICATES;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_INTRONIC;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_READ_THROUGH;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_TOTAL;
import static com.hartwig.hmftools.svtools.rna_expression.GeneReadData.GC_TRANS_SUPPORTING;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.markRegionBases;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.validRegionMatchType;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.validTranscriptType;
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

    // state relating to the current gene
    private GeneReadData mCurrentGene;
    private final List<String> mDiscardedReads;

    private int mGeneReadCount;
    private int mTotalBamReadCount;

    private static final int DEFAULT_MIN_MAPPING_QUALITY = 1;

    private final Map<String,ReadRecord> mFragmentReads;
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

        mWriter = null;
    }

    public boolean validReader() { return mSamReader != null; }

    public void close()
    {
        LOGGER.info("read {} total BAM records", mTotalBamReadCount);

        closeBufferedWriter(mWriter);
    }

    public void readBamCounts(final GeneReadData geneReadData, final GenomeRegion genomeRegion)
    {
        mFragmentReads.clear();
        mDiscardedReads.clear();

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
            mCurrentGene.addCount(GC_TOTAL, mGeneReadCount / 2);
        }
    }

    public void readBamCounts(final GenomeRegion genomeRegion, final Consumer<SAMRecord> consumer)
    {
        SAMSlicer samSlicer = new SAMSlicer(DEFAULT_MIN_MAPPING_QUALITY, Lists.newArrayList(genomeRegion));
        samSlicer.slice(mSamReader, consumer);
    }

    private void processSamRecord(@NotNull final SAMRecord record)
    {
        if(record.getDuplicateReadFlag())
        {
            if(record.getFirstOfPairFlag())
                ++mCurrentGene.getCounts()[GC_DUPLICATES];

            if(!mConfig.KeepDuplicates)
                return;
        }

        ++mTotalBamReadCount;
        ++mGeneReadCount;

        if(mConfig.GeneStatsOnly)
        {
            recordFragmentLength(record);
            return;
        }

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
            mCurrentGene.addCount(GC_TOTAL, 1);
            mCurrentGene.addCount(GC_CHIMERIC, 1);
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
                -
        */

        boolean r1OutsideGene = read1.PosStart > mCurrentGene.GeneData.GeneEnd || read1.PosEnd < mCurrentGene.GeneData.GeneStart;
        boolean r2OutsideGene = read2.PosStart > mCurrentGene.GeneData.GeneEnd || read2.PosEnd < mCurrentGene.GeneData.GeneStart;

        if(r1OutsideGene && r2OutsideGene)
            return;

        mCurrentGene.addCount(GC_TOTAL, 1);

        if(read1.localInversion() || read2.localInversion())
        {
            mCurrentGene.addCount(GC_CHIMERIC, 1);
            return;
        }

        if(read1.getMappedRegions().isEmpty() && read2.getMappedRegions().isEmpty())
            return;

        final Map<String,TransMatchType> firstReadTransTypes = read1.getTranscriptClassifications();

        final Map<String,TransMatchType> secondReadTransTypes = read2.getTranscriptClassifications();

        // first find valid transcripts in both reads
        final List<String> firstReadValidTrans = firstReadTransTypes.entrySet().stream()
                .filter(x -> validTranscriptType(x.getValue()))
                .map(x -> x.getKey()).collect(Collectors.toList());

        final List<String> validTranscripts = secondReadTransTypes.entrySet().stream()
                .filter(x -> validTranscriptType(x.getValue()))
                .filter(x -> firstReadValidTrans.contains(x.getKey()))
                .map(x -> x.getKey()).collect(Collectors.toList());

        boolean isLongFragment = read1.fragmentInsertSize() > mConfig.LongFragmentLimit;

        // now mark all other transcripts which aren't valid either due to the read pair
        if(!validTranscripts.isEmpty())
        {
            firstReadTransTypes.entrySet().stream()
                    .filter(x -> validTranscriptType(x.getValue()))
                    .filter(x -> !validTranscripts.contains(x.getKey()))
                    .forEach(x -> x.setValue(OTHER_TRANS));

            secondReadTransTypes.entrySet().stream()
                    .filter(x -> validTranscriptType(x.getValue()))
                    .filter(x -> !validTranscripts.contains(x.getKey()))
                    .forEach(x -> x.setValue(OTHER_TRANS));
        }
        else
        {
            firstReadTransTypes.entrySet().stream()
                    .filter(x -> validTranscriptType(x.getValue()))
                    .forEach(x -> x.setValue(ALT));

            secondReadTransTypes.entrySet().stream()
                    .filter(x -> validTranscriptType(x.getValue()))
                    .forEach(x -> x.setValue(ALT));


            if(r1OutsideGene || r2OutsideGene)
                mCurrentGene.addCount(GC_READ_THROUGH, 1);
            else
                mCurrentGene.addCount(GC_ALT, 1);
        }

        // finally record valid read info against each region now that it is known
        if(!validTranscripts.isEmpty())
        {
            mCurrentGene.addCount(GC_TRANS_SUPPORTING, 1);

            // now record the bases covered by the read in these matched regions
            final List<long[]> commonMappings = deriveCommonRegions(read1.getMappedRegionCoords(), read2.getMappedRegionCoords());

            final List<RegionReadData> regions = read1.getMappedRegions().entrySet().stream()
                    .filter(x -> validRegionMatchType(x.getValue()))
                    .map(x -> x.getKey()).collect(Collectors.toList());

            final List<RegionReadData> regions2 = read2.getMappedRegions().entrySet().stream()
                    .filter(x -> validRegionMatchType(x.getValue()))
                    .map(x -> x.getKey()).collect(Collectors.toList());

            for(RegionReadData region : regions2)
            {
                if (!regions.contains(region))
                    regions.add(region);
            }

            regions.forEach(x -> markRegionBases(commonMappings, x));

            boolean isUniqueTrans = validTranscripts.size() == 1;

            for (final String trans : validTranscripts)
            {
                mCurrentGene.addTranscriptReadMatch(trans, isUniqueTrans);

                // keep track of which regions have been allocated from this fragment as a whole, so not counting each read separately
                final List<RegionReadData> processedRegions = Lists.newArrayList();

                processValidTranscript(trans, read1, processedRegions, isUniqueTrans);
                processValidTranscript(trans, read2, processedRegions, isUniqueTrans);
            }
        }

        if(mConfig.WriteReadData)
        {
            writeReadData(0, read1);
            writeReadData(1, read2);
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
            mCurrentGene.addCount(GC_INTRONIC, 1);
            mCurrentGene.addCount(GC_TOTAL, 1);
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

        recordFragmentLength(read.samRecord());
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

    private void recordFragmentLength(SAMRecord record)
    {
        if(!mConfig.WriteFragmentLengths || record == null)
            return;

        if(!record.getMateReferenceName().equals(record.getReferenceName())
        || record.getMateNegativeStrandFlag() != record.getReadNegativeStrandFlag())
            return;

        int fragmentSize = record.getInferredInsertSize();
        if (fragmentSize > 0)
        {
            mCurrentGene.addFragmentLength(fragmentSize);
        }
    }

    private void writeReadData(int readIndex, final ReadRecord read)
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
                mWriter.write(",TransId,TransClass,ExonRank,ExonStart,ExonEnd,MatchType");
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

                    mWriter.write(String.format(",%s,%s,%d,%d,%d,%s",
                            trans, transType, region.getExonRank(trans), region.start(), region.end(), matchType));

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
