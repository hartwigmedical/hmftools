package com.hartwig.hmftools.svtools.rna_expression;

import static java.lang.Math.min;

import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.io.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.validRegionMatchType;
import static com.hartwig.hmftools.svtools.rna_expression.ReadRecord.validTranscriptType;
import static com.hartwig.hmftools.svtools.rna_expression.RegionMatchType.INTRONIC;
import static com.hartwig.hmftools.svtools.rna_expression.RnaExpConfig.MAX_READ_COUNT;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.ALT;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.OTHER_TRANS;
import static com.hartwig.hmftools.svtools.rna_expression.TransMatchType.SPLICE_JUNCTION;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.hotspot.SAMSlicer;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;

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
    private final List<ReadRecord> mReadRecords;
    private int mBamRecordCount;
    private GeneReadData mCurrentGene;

    private static final int DEFAULT_MIN_MAPPING_QUALITY = 1;
    private static final double MIN_BASE_MATCH_PERC = 0.9;

    private final Map<String,ReadRecord> mFragmentReads;
    private BufferedWriter mWriter;

    private static final Logger LOGGER = LogManager.getLogger(RnaBamReader.class);

    public RnaBamReader(final RnaExpConfig config)
    {
        mConfig = config;

        mReadRecords = Lists.newArrayList();
        mBamRecordCount = 0;
        mCurrentGene = null;
        mFragmentReads = Maps.newHashMap();

        mSamReader = SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile));
        mWriter = null;
    }

    public void close() { closeBufferedWriter(mWriter); }

    public void readBamCounts(final GeneReadData geneReadData, final GenomeRegion genomeRegion)
    {
        mReadRecords.clear();
        mFragmentReads.clear();
        mBamRecordCount = 0;

        mCurrentGene = geneReadData;

        SAMSlicer samSlicer = new SAMSlicer(DEFAULT_MIN_MAPPING_QUALITY, Lists.newArrayList(genomeRegion));
        samSlicer.slice(mSamReader, this::processSamRecord);
    }

    public void readBamCounts(final GenomeRegion genomeRegion, final Consumer<SAMRecord> consumer)
    {
        SAMSlicer samSlicer = new SAMSlicer(DEFAULT_MIN_MAPPING_QUALITY, Lists.newArrayList(genomeRegion));
        samSlicer.slice(mSamReader, consumer);
    }

    private void processSamRecord(@NotNull final SAMRecord record)
    {
        // check for records which overlap an exonic regions
        boolean exonOverlap = mCurrentGene.getExonRegions().stream()
                .anyMatch(x -> !(record.getEnd() < x.start() || record.getStart() > x.Region.end()));

        if(!exonOverlap)
        {
            checkIntronicRegions(record);
            return;
        }

        ++mBamRecordCount;

        if(mConfig.ReadCountLimit > 0 && mReadRecords.size() >= mConfig.ReadCountLimit)
            return;

        mReadRecords.add(ReadRecord.from(record));
    }

    public void analyseReads()
    {
        mCurrentGene.setTotalReadCount(mBamRecordCount);

        if(mBamRecordCount >= MAX_READ_COUNT)
        {
            LOGGER.warn("gene({}) readCount({}) exceeds max read count", mCurrentGene.GeneData.GeneName, mBamRecordCount);
            // return;
        }

        // cache reference bases for comparison with read bases
        for(RegionReadData region : mCurrentGene.getExonRegions())
        {
            final String regionRefBases = mConfig.RefFastaSeqFile.getSubsequenceAt(
                    region.chromosome(), region.start(), region.end()).getBaseString();

            region.setRefBases(regionRefBases);
        }

        // for each record find all exons with an overlap
        // skip records if either end isn't in one of the exons for this gene

        for(final ReadRecord read : mReadRecords)
        {
            // the read is fully within the exon
            List<RegionReadData> overlappingRegions = mCurrentGene.getExonRegions().stream()
                    .filter(x -> read.overlapsMappedReads(x.Region.start(), x.Region.end()))
                    .collect(Collectors.toList());

            if(overlappingRegions.isEmpty())
                continue;

            // look at all matched reads within the context of a transcript
            read.processOverlappingRegions(overlappingRegions);

            ReadRecord otherRead = checkFragmentRead(read);

            if(otherRead != null)
                processFragmentReads(read, otherRead);
        }

        if(!mFragmentReads.isEmpty())
        {
            LOGGER.debug("gene({}) has {} unmatched reads", mCurrentGene.GeneData.GeneName, mFragmentReads.size());
            mFragmentReads.clear();
        }
    }

    public void processFragmentReads(@NotNull final ReadRecord read1, @NotNull final ReadRecord read2)
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
        }

        // finally record valid read info against each region now that it is known
        for(final String trans : validTranscripts)
        {
            mCurrentGene.addTranscriptReadMatch(trans);

            for(int i = 0; i <= 1; ++i)
            {
                ReadRecord read = (i == 0) ? read1 : read2;

                List<RegionReadData> regions = read.getMappedRegions().entrySet().stream()
                        .filter(x -> x.getKey().hasTransId(trans))
                        .filter(x -> validRegionMatchType(x.getValue()))
                        .map(x -> x.getKey()).collect(Collectors.toList());

                regions.forEach(x -> x.addTranscriptReadMatch(trans));

                // now record the bases covered by the read in these matched regions

                // any adjacent reads can record a splice junction count

                if(regions.size() > 1 && read.getTranscriptClassification(trans) == SPLICE_JUNCTION)
                {
                    for(int r1 = 0; r1 < regions.size() - 1; ++r1)
                    {
                        RegionReadData region1 = regions.get(r1);

                        for(int r2 = r1 + 1; r2 < regions.size(); ++r2)
                        {
                            RegionReadData region2 = regions.get(r2);

                            if(region1.getPostRegions().contains(region2))
                            {
                                region1.addTranscriptJunctionMatch(trans, SE_END);
                                region2.addTranscriptJunctionMatch(trans, SE_START);
                            }
                            else if(region1.getPreRegions().contains(region2))
                            {
                                region1.addTranscriptJunctionMatch(trans, SE_START);
                                region2.addTranscriptJunctionMatch(trans, SE_END);
                            }
                        }
                    }
                }
            }
        }

        if(mConfig.WriteReadData)
        {
            writeReadData(read1);
            writeReadData(read2);
        }
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

    private void checkIntronicRegions(final SAMRecord record)
    {
        if(record.getCigar() == null)
            return;

        if(record.getCigar().containsOperator(CigarOperator.N) || !record.getCigar().containsOperator(CigarOperator.M))
            return;

        RegionReadData intronReadData = mCurrentGene.getIntronRegions().stream()
                .filter(x -> record.getStart() >= x.Region.start() && record.getStart() <= x.Region.end())
                .findFirst().orElse(null);

        if(intronReadData != null)
        {
            if(mConfig.AllTranscripts && intronReadData.getRefRegions().size() == 1)
            {
                // only record intronic reads if they are unique to a transcript
                intronReadData.addMatchedRead(INTRONIC);
            }

            if(record.getInferredInsertSize() > 0)
            {
                // cache reads likely to map into an exon
                if(record.getStart() - record.getInferredInsertSize() <= intronReadData.start()
                || record.getEnd() + record.getInferredInsertSize() >= intronReadData.start())
                {
                    ReadRecord read = ReadRecord.from(record);
                    ReadRecord otherRead = checkFragmentRead(read);

                    if(otherRead != null)
                        processFragmentReads(read, otherRead);
                }
            }
        }

        if(mConfig.WriteFragmentLengths)
        {
            int fragmentSize = record.getInferredInsertSize();
            if (fragmentSize > 0 && record.getMateReferenceName().equals(record.getReferenceName()))
            {
                mCurrentGene.addFragmentLength(fragmentSize);
            }
        }
    }

    private ReadRecord checkFragmentRead(ReadRecord read)
    {
        if(read.samRecord() != null)
        {
            if(!read.samRecord().getMateReferenceName().equals(read.Chromosome)
            || read.samRecord().getMateReferenceIndex() == null)
            {
                return null;
            }
        }

        ReadRecord otherRead = mFragmentReads.get(read.Id);

        if(otherRead != null)
            return otherRead;

        mFragmentReads.put(read.Id, read);
        return null;
    }

    private void writeReadData(final ReadRecord read)
    {
        if(mConfig.OutputDir.isEmpty())
            return;

        try
        {
            if(mWriter == null)
            {
                final String outputFileName = mConfig.OutputDir + "RNA_READ_DATA.csv";

                mWriter = createBufferedWriter(outputFileName, false);
                mWriter.write("GeneId,GeneName,ReadId,Chromosome,PosStart,PosEnd,Cigar");
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

                    mWriter.write(String.format("%s,%s,%s,%s,%d,%d,%s",
                            mCurrentGene.GeneData.GeneId, mCurrentGene.GeneData.GeneName, read.Id,
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

    public static int findStringOverlaps(final String str1, final String str2)
    {
        if(str1.length() == 0 || str2.length() == 0)
            return 0;

        int matched = 0;
        int i = 0;
        int j = 0;
        int mismatchIndex = -1;

        // first compare bases at same indices, making note of the first difference if there is one
        while(i < str1.length() && j < str2.length())
        {
            if (str1.charAt(i) == str2.charAt(j))
                ++matched;
            else if(mismatchIndex == -1)
                mismatchIndex = i;

            ++i;
            ++j;
        }

        if(matched > MIN_BASE_MATCH_PERC * min(str1.length(), str2.length()))
            return matched;

        i = j = mismatchIndex;
        matched = mismatchIndex;

        while(i < str1.length() && j < str2.length())
        {
            if(str1.charAt(i) == str2.charAt(j))
            {
                ++i;
                ++j;
                ++matched;
                continue;
            }

            // search ahead in each string in turn for the next short matching sequence
            int startI = i;
            boolean seqFound = false;
            for(; i < str1.length() - 2 && j < str2.length() - 2; ++i)
            {
                if(str1.charAt(i) == str2.charAt(j) && str1.charAt(i+1) == str2.charAt(j+1) && str1.charAt(i+2) == str2.charAt(j+2))
                {
                    seqFound = true;
                    break;
                }
            }

            if(seqFound)
                continue;

            i = startI;

            for(; i < str1.length() - 2 && j < str2.length() - 2; ++j)
            {
                if(str1.charAt(i) == str2.charAt(j) && str1.charAt(i+1) == str2.charAt(j+1) && str1.charAt(i+2) == str2.charAt(j+2))
                {
                    seqFound = true;
                    break;
                }
            }

            if(!seqFound)
                break;
        }

        return matched;
    }

    @VisibleForTesting
    public void addReadRecords(final GeneReadData geneReadData, final List<ReadRecord> readRecords)
    {
        mCurrentGene = geneReadData;
        mBamRecordCount += readRecords.size();
        mReadRecords.clear();;
        mReadRecords.addAll(readRecords);
    }


}
