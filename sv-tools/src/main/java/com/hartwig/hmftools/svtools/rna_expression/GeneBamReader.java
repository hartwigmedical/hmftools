package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.ALT;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.CHIMERIC;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.DUPLICATE;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.READ_THROUGH;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.TOTAL;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.TRANS_SUPPORTING;
import static com.hartwig.hmftools.svtools.rna_expression.GeneMatchType.UNSPLICED;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.calcFragmentLength;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.getUniqueValidRegion;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.hasSkippedExons;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.markRegionBases;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.validRegionMatchType;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.validTranscriptType;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.EXON_INTRON;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.deriveCommonRegions;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpUtils.positionsOverlap;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.OTHER_TRANS;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.SPLICE_JUNCTION;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class GeneBamReader
{
    private final RnaExpConfig mConfig;
    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    // state relating to the current gene
    private GeneReadData mCurrentGene;
    private final Map<String,ReadRecord> mFragmentReads; // delay processing of read until both have been read

    private int mGeneReadCount;
    private int mTotalBamReadCount;

    public static final int DEFAULT_MIN_MAPPING_QUALITY = 1;

    private final Map<Long,List<long[]>> mDuplicateCache;
    private final List<String> mDuplicateReadIds;

    private final List<TranscriptComboData> mTransComboData;
    private final AltSpliceJunctionFinder mAltSpliceJunctionFinder;
    private final FragmentSizeCalcs mFragmentSizeCalc;

    private final BufferedWriter mReadDataWriter;

    private static final Logger LOGGER = LogManager.getLogger(GeneBamReader.class);

    public GeneBamReader(final RnaExpConfig config, final ResultsWriter resultsWriter, final FragmentSizeCalcs fragmentSizeCalc)
    {
        mConfig = config;

        mCurrentGene = null;
        mFragmentReads = Maps.newHashMap();
        mTransComboData = Lists.newArrayList();

        mGeneReadCount = 0;
        mTotalBamReadCount = 0;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        mBamSlicer = new BamSlicer(DEFAULT_MIN_MAPPING_QUALITY, false);

        mDuplicateCache = Maps.newHashMap();
        mDuplicateReadIds = Lists.newArrayList();

        mReadDataWriter = resultsWriter != null ? resultsWriter.getReadDataWriter() : null;

        mAltSpliceJunctionFinder = new AltSpliceJunctionFinder(
                mConfig, mSamReader, resultsWriter != null ? resultsWriter.getAltSpliceJunctionWriter() : null);

        mFragmentSizeCalc = fragmentSizeCalc;
    }

    public static GeneBamReader from(final RnaExpConfig config)
    {
        return new GeneBamReader(config, null, null);
    }

    public int totalBamCount() { return mTotalBamReadCount; }

    public void readBamCounts(final GeneReadData geneReadData, final GenomeRegion genomeRegion)
    {
        mFragmentReads.clear();

        clearDuplicates();

        mCurrentGene = geneReadData;
        mGeneReadCount = 0;
        mTransComboData.clear();
        mAltSpliceJunctionFinder.setGeneData(geneReadData);

        mBamSlicer.slice(mSamReader, Lists.newArrayList(genomeRegion), this::processSamRecord);

        LOGGER.debug("gene({}) bamReadCount({})", mCurrentGene.GeneData.GeneName, mGeneReadCount);

        if(mConfig.GeneStatsOnly)
        {
            mCurrentGene.addCount(TOTAL, mGeneReadCount / 2);
        }

        mAltSpliceJunctionFinder.recordDepthCounts();
        mAltSpliceJunctionFinder.writeAltSpliceJunctions();
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

        if(read.isTranslocation())
        {
            if(read.samRecord().getFirstOfPairFlag())
            {
                mCurrentGene.addCount(TOTAL, 1);
                mCurrentGene.addCount(CHIMERIC, 1);
            }
            return;
        }

        boolean exonOverlap = mCurrentGene.getExonRegions().stream()
                .anyMatch(x -> positionsOverlap(read.PosStart, read.PosEnd, x.start(), x.end()));

        /*
        if(!exonOverlap)
        {
            boolean outsideGene = !positionsOverlap(read.PosStart, read.PosEnd, mCurrentGene.GeneData.GeneStart, mCurrentGene.GeneData.GeneEnd);

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
        */

        if(exonOverlap)
        {
            List<RegionReadData> overlappingRegions = mCurrentGene.findOverlappingRegions(read);

            if (!overlappingRegions.isEmpty())
            {
                // look at all matched reads within the context of a transcript
                read.processOverlappingRegions(overlappingRegions);
            }
        }

        checkFragmentRead(read);
    }

    private boolean checkFragmentRead(ReadRecord read)
    {
        // check if the 2 reads from a fragment exist and if so handle them a pair, returning true
        if(read.samRecord() != null)
        {
            if(read.isTranslocation() || read.samRecord().getMateReferenceIndex() == null)
                return false;
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

        boolean r1OutsideGene = read1.PosStart > mCurrentGene.GeneData.GeneEnd || read1.PosEnd < mCurrentGene.GeneData.GeneStart;
        boolean r2OutsideGene = read2.PosStart > mCurrentGene.GeneData.GeneEnd || read2.PosEnd < mCurrentGene.GeneData.GeneStart;

        if(r1OutsideGene && r2OutsideGene)
            return;

        mCurrentGene.addCount(TOTAL, 1);

        if(read1.isLocalInversion() || read2.isLocalInversion())
        {
            mCurrentGene.addCount(CHIMERIC, 1);
            return;
        }

        if(read1.getMappedRegions().isEmpty() && read2.getMappedRegions().isEmpty())
        {
            // fully intronic read
            processIntronicReads(read1, read2);
            return;
        }

        final Map<Integer,TransMatchType> firstReadTransTypes = read1.getTranscriptClassifications();

        final Map<Integer,TransMatchType> secondReadTransTypes = read2.getTranscriptClassifications();

        // first find valid transcripts in both reads
        final List<Integer> validTranscripts = Lists.newArrayList();
        final List<Integer> invalidTranscripts = Lists.newArrayList();
        int calcFragmentLength = calcFragmentLength(read1, read2);
        boolean validFragmentLength = calcFragmentLength <= mConfig.MaxFragmentLength;

        final List<RegionReadData> validRegions = getUniqueValidRegion(read1, read2);

        if(mConfig.RunValidations)
        {
            for(RegionReadData region : validRegions)
            {
                if(validRegions.stream().filter(x -> x == region).count() > 1)
                {
                    LOGGER.error("repeated exon region");
                }
            }
        }

        for(Map.Entry<Integer,TransMatchType> entry : firstReadTransTypes.entrySet())
        {
            int transId = entry.getKey();

            if(validFragmentLength && validTranscriptType(entry.getValue()))
            {
                if(secondReadTransTypes.containsKey(transId) && validTranscriptType(secondReadTransTypes.get(transId)))
                {
                    if(!hasSkippedExons(validRegions, transId, mConfig.MaxFragmentLength))
                    {
                        validTranscripts.add(transId);
                        continue;
                    }
                }
            }

            if(!invalidTranscripts.contains(transId))
                invalidTranscripts.add(transId);
        }

        for(Integer transId : secondReadTransTypes.keySet())
        {
            if(!validTranscripts.contains(transId) && !invalidTranscripts.contains(transId))
                invalidTranscripts.add(transId);
        }

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
                mAltSpliceJunctionFinder.evaluateFragmentReads(read1, read2, invalidTranscripts);
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

            if(mConfig.RunValidations)
            {
                for(long[] readRegion : commonMappings)
                {
                    if(commonMappings.stream().filter(x -> x[SE_START] == readRegion[SE_START] && x[SE_END] == readRegion[SE_END]).count() > 1)
                    {
                        LOGGER.error("repeated read region");
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
                    transMatchType = FragmentMatchType.LONG;

                    if(comboTransMatchType != FragmentMatchType.SPLICED)
                        comboTransMatchType = FragmentMatchType.LONG;
                }
                else
                {
                    transMatchType = FragmentMatchType.SHORT;
                }

                mCurrentGene.addTranscriptReadMatch(transId, isUniqueTrans, transMatchType);

                // keep track of which regions have been allocated from this fragment as a whole, so not counting each read separately
                final List<RegionReadData> processedRegions = Lists.newArrayList();

                processValidTranscript(transId, read1, processedRegions, isUniqueTrans);
                processValidTranscript(transId, read2, processedRegions, isUniqueTrans);
            }

            addTransComboData(validTranscripts, comboTransMatchType);
        }

        mCurrentGene.addCount(geneReadType, 1);

        if(mConfig.WriteReadData && mReadDataWriter != null)
        {
            writeReadData(mReadDataWriter, mCurrentGene, 0, read1, geneReadType, validTranscripts.size(), calcFragmentLength);
            writeReadData(mReadDataWriter, mCurrentGene, 1, read2, geneReadType, validTranscripts.size(), calcFragmentLength);
        }
    }

    public List<TranscriptComboData> getTransComboData() { return mTransComboData; }

    private void addTransComboData(final List<Integer> transcripts, FragmentMatchType transMatchType)
    {
        TranscriptComboData transComboCounts = mTransComboData.stream()
                .filter(x -> x.matches(transcripts)).findFirst().orElse(null);

        if(transComboCounts == null)
        {
            transComboCounts = new TranscriptComboData(transcripts);
            mTransComboData.add(transComboCounts);
        }

        transComboCounts.addCounts(transMatchType, 1);
    }

    private void processValidTranscript(
            int transId, final ReadRecord read, final List<RegionReadData> processedRegions, boolean isUniqueTrans)
    {
        List<RegionReadData> regions = read.getMappedRegions().entrySet().stream()
                .filter(x -> x.getKey().hasTransId(transId))
                .filter(x -> validRegionMatchType(x.getValue()))
                .map(x -> x.getKey()).collect(Collectors.toList());

        for(RegionReadData region : regions)
        {
            if (!processedRegions.contains(region))
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

    private void processIntronicReads(final ReadRecord read1, final ReadRecord read2)
    {
        if(read1.Cigar.containsOperator(CigarOperator.N) || read2.Cigar.containsOperator(CigarOperator.N))
        {
            if(mCurrentGene.overlapsOtherGeneExon(read1.PosStart, read1.PosEnd) || mCurrentGene.overlapsOtherGeneExon(read2.PosStart, read2.PosEnd))
                return;

            mCurrentGene.addCount(ALT, 1);
            mAltSpliceJunctionFinder.evaluateFragmentReads(read1, read2, Lists.newArrayList());
            return;
        }

        long fragMinPos = min(read1.PosStart, read2.PosStart);
        long fragMaxPos = max(read1.PosEnd, read2.PosEnd);

        if(mCurrentGene.overlapsOtherGeneExon(fragMinPos, fragMaxPos))
            return;

        mCurrentGene.addCount(UNSPLICED, 1);

        // both reads will report the same fragment length (ie insert size for non-split reads)
        if(mFragmentSizeCalc != null)
            mFragmentSizeCalc.addFragmentLength(read1.samRecord());
    }

    /*
    private void checkIntronicRegions(final ReadRecord read)
    {
        if(read.Cigar == null)
            return;

        if(!read.Cigar.containsOperator(CigarOperator.M) || read.isLocalInversion())
            return;

        if(read.Cigar.containsOperator(CigarOperator.N))
        {
            if(!mCurrentGene.overlapsOtherGeneExon(read.PosStart, read.PosEnd))
            {
                mAltSpliceJunctionFinder.evaluateIntronicRead(read);
            }

            return;
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

        // only count this read as intronic if it doesn't overlap with other genes' exons
        long fragMinPos = min(otherReadStartPos, read.PosStart);
        long fragMaxPos = max(otherReadEndPos, read.PosEnd);

        if(mCurrentGene.overlapsOtherGeneExon(fragMinPos, fragMaxPos))
            return;

        if(read.PosStart < otherReadStartPos)
        {
            mCurrentGene.addCount(UNSPLICED, 1);
            mCurrentGene.addCount(TOTAL, 1);
        }
    }
    */

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

    public static BufferedWriter createReadDataWriter(final RnaExpConfig config)
    {
        try
        {
            final String outputFileName = config.formOutputFile("read_data.csv");

            BufferedWriter writer = createBufferedWriter(outputFileName, false);
            writer.write("GeneId,GeneName,ReadIndex,ReadId,Chromosome,PosStart,PosEnd,Cigar,InsertSize,FragLength");
            writer.write(",GeneClass,TransId,TransClass,ValidTrans,ExonRank,ExonStart,ExonEnd,RegionClass");
            writer.newLine();
            return writer;
        }
        catch (IOException e)
        {
            LOGGER.error("failed to create read data writer: {}", e.toString());
            return null;
        }
    }

    private synchronized static void writeReadData(
            final BufferedWriter writer, final GeneReadData geneReadData, int readIndex, final ReadRecord read,
            GeneMatchType geneReadType, int validTranscripts, int calcFragmentLength)
    {
        if(read.getTranscriptClassifications().isEmpty())
            return;

        try
        {
            for(Map.Entry<Integer,TransMatchType> entry : read.getTranscriptClassifications().entrySet())
            {
                int transId = entry.getKey();
                TransMatchType transType = entry.getValue();

                for(Map.Entry<RegionReadData,RegionMatchType> rEntry : read.getMappedRegions().entrySet())
                {
                    RegionReadData region = rEntry.getKey();
                    RegionMatchType matchType = rEntry.getValue();

                    if(!region.hasTransId(transId))
                        continue;

                    writer.write(String.format("%s,%s,%d,%s",
                            geneReadData.GeneData.GeneId, geneReadData.GeneData.GeneName, readIndex, read.Id));

                    writer.write(String.format(",%s,%d,%d,%s,%d,%d",
                            read.Chromosome, read.PosStart, read.PosEnd, read.Cigar.toString(),
                            read.fragmentInsertSize(), calcFragmentLength));

                    writer.write(String.format(",%s,%d,%s,%s,%d,%d,%d,%s",
                            geneReadType, transId, transType, validTranscripts,
                            region.getExonRank(transId), region.start(), region.end(), matchType));

                    writer.newLine();
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
