package com.hartwig.hmftools.fastqtools.biomodalcollapse;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.fastqtools.FastqCommon.FQ_LOGGER;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalCollapseUtil.sanatizeQualString;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalConstants.COMPARE_DELIMITER;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalConstants.COMPARE_HEADERS;
import static com.hartwig.hmftools.fastqtools.biomodalcollapse.BiomodalConstants.MISSING_BASE;
import static com.hartwig.hmftools.common.bam.CigarUtils.leftSoftClipLength;
import static com.hartwig.hmftools.common.bam.CigarUtils.rightSoftClipLength;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.SUPP_POS_STRAND;
import static com.hartwig.hmftools.common.bam.SupplementaryReadData.extractAlignments;
import static com.hartwig.hmftools.common.utils.PerformanceCounter.runTimeMinsStr;

import static htsjdk.samtools.CigarOperator.D;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.M;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.SortedMap;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import com.google.common.collect.Lists;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BiomodalCompare
{
    private static final int SHORT_INV_DISTANCE_CUTOFF = 100;

    private static final int LOW_QUAL_CUTOFF = 19;
    private static final int MID_QUAL_CUTOFF = 30;

    private static class ReadPair
    {
        public SAMRecord RefRead;
        public SAMRecord NewRead;

        public ReadPair()
        {
            RefRead = null;
            NewRead = null;
        }
    }

    private final BiomodalCompareConfig mConfig;
    private final SortedMap<String, ReadPair> mReadPairMap;

    public BiomodalCompare(final BiomodalCompareConfig config)
    {
        mConfig = config;
        mReadPairMap = Maps.newTreeMap();
    }

    public void loadReads(final SamReader samReader, boolean isRef)
    {
        SAMRecordIterator iter = samReader.iterator();
        while(iter.hasNext())
        {
            SAMRecord read = iter.next();
            if(read.isSecondaryOrSupplementary())
            {
                continue;
            }

            String readName = read.getReadName();
            ReadPair readPair = mReadPairMap.get(readName);
            if(readPair == null)
            {
                readPair = new ReadPair();
                mReadPairMap.put(readName, readPair);
            }

            if(isRef)
            {
                readPair.RefRead = read;
            }
            else
            {
                readPair.NewRead = read;
            }
        }
    }

    public void run()
    {
        FQ_LOGGER.info("starting BiomodalCompare...");
        long startTimeMs = System.currentTimeMillis();

        SamReaderFactory readerFactory = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);

        try(SamReader refReader = readerFactory.open(new File(mConfig.RefBam));
                SamReader newReader = readerFactory.open(new File(mConfig.NewBam));
                BufferedWriter compareWriter = new BufferedWriter(new FileWriter(mConfig.OutputTsv)))
        {
            compareWriter.write(Arrays.stream(COMPARE_HEADERS).collect(Collectors.joining(COMPARE_DELIMITER)));
            compareWriter.newLine();

            loadReads(refReader, true);
            loadReads(newReader, false);

            for(ReadPair readPair : mReadPairMap.values())
            {
                processReadPair(compareWriter, readPair.RefRead, readPair.NewRead);
            }
        }
        catch(FileNotFoundException e)
        {
            throw new RuntimeException(e);
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }

        FQ_LOGGER.info("BiomodalCompare complete, mins({})", runTimeMinsStr(startTimeMs));
    }

    private static class AlignmentInfo
    {
        public final String Chromosome;
        public final int Pos;
        public final int AlignmentLength;
        public final int MatchCount;
        public final int LeftSoftClipLength;
        public final int RightSoftClipLength;
        public final int InsertCount;
        public final int DeleteCount;

        public AlignmentInfo(final String chromosome, int pos, int alignmentLength, int matchCount, int leftSoftClipLength,
                int rightSoftClipLength, int insertCount, int deleteCount)
        {
            Chromosome = chromosome;
            Pos = pos;
            AlignmentLength = alignmentLength;
            MatchCount = matchCount;
            LeftSoftClipLength = leftSoftClipLength;
            RightSoftClipLength = rightSoftClipLength;
            InsertCount = insertCount;
            DeleteCount = deleteCount;
        }
    }

    @Nullable
    private static AlignmentInfo getAlignmentInfo(@Nullable final SAMRecord read)
    {
        if(read == null)
        {
            return null;
        }

        if(read.getReadUnmappedFlag())
        {
            return null;
        }

        int matchCount = 0;
        int insertCount = 0;
        int deleteCount = 0;
        for(CigarElement el : read.getCigar().getCigarElements())
        {
            if(el.getOperator() == M)
            {
                matchCount += el.getLength();
            }
            else if(el.getOperator() == I)
            {
                insertCount += el.getLength();
            }
            else if(el.getOperator() == D)
            {
                deleteCount += el.getLength();
            }
        }

        return new AlignmentInfo(
                read.getContig(),
                read.getAlignmentStart(),
                read.getAlignmentEnd() - read.getAlignmentStart() + 1,
                matchCount,
                leftSoftClipLength(read),
                rightSoftClipLength(read),
                insertCount,
                deleteCount
        );
    }

    private static boolean isShortInvSupp(final SAMRecord read, final SupplementaryReadData suppData)
    {
        if(read.getReadUnmappedFlag())
        {
            return false;
        }

        if(!read.getContig().equals(suppData.Chromosome))
        {
            return false;
        }

        if((suppData.Strand != SUPP_POS_STRAND) == read.getReadNegativeStrandFlag())
        {
            return false;
        }

        int posDiff = abs(read.getAlignmentStart() - suppData.Position);
        return posDiff < SHORT_INV_DISTANCE_CUTOFF;
    }

    private static class SuppInfo
    {
        public final List<SupplementaryReadData> ShortInvSupps;
        public final List<SupplementaryReadData> OtherSupps;

        public SuppInfo()
        {
            ShortInvSupps = Lists.newArrayList();
            OtherSupps = Lists.newArrayList();
        }
    }

    @Nullable
    private static SuppInfo getSuppInfo(@Nullable final SAMRecord read)
    {
        if(read == null || read.getReadUnmappedFlag())
        {
            return null;
        }

        List<SupplementaryReadData> supps = extractAlignments(read);
        if(supps == null)
        {
            return null;
        }

        SuppInfo suppInfo = new SuppInfo();
        for(SupplementaryReadData supp : supps)
        {
            if(isShortInvSupp(read, supp))
            {
                suppInfo.ShortInvSupps.add(supp);
            }
            else
            {
                suppInfo.OtherSupps.add(supp);
            }
        }

        return suppInfo;
    }

    private static class MDOperator
    {
        public final CigarOperator Operator;
        public final Character RefBase;

        public MDOperator(final CigarOperator operator, @Nullable final Character refBase)
        {
            Operator = operator;
            RefBase = refBase;
        }
    }

    private static String parseIntPrefix(final String s, int startIndex)
    {
        StringBuilder intPrefix = new StringBuilder();
        for(int i = startIndex; i < s.length(); i++)
        {
            char c = s.charAt(i);
            if(c >= '0' && c <= '9')
            {
                intPrefix.append(c);
                continue;
            }

            break;
        }

        return intPrefix.toString();
    }

    private static String parseNonIntPrefix(final String s, int startIndex)
    {
        StringBuilder prefix = new StringBuilder();
        for(int i = startIndex; i < s.length(); i++)
        {
            char c = s.charAt(i);
            if(c >= '0' && c <= '9')
            {
                break;
            }

            prefix.append(c);
        }

        return prefix.toString();
    }

    @Nullable
    private static List<MDOperator> parseMDTag(@Nullable final SAMRecord read)
    {
        if(read == null || read.getReadUnmappedFlag())
        {
            return null;
        }

        String md = read.getStringAttribute("MD");
        if(md == null)
        {
            return null;
        }

        List<MDOperator> ops = Lists.newArrayList();
        int i = 0;
        while(i < md.length())
        {
            String exactMatchCountStr = parseIntPrefix(md, i);
            i += exactMatchCountStr.length();
            int exactMatchCount = Integer.parseInt(exactMatchCountStr);
            for(int j = 0; j < exactMatchCount; j++)
            {
                ops.add(new MDOperator(M, null));
            }

            if(i >= md.length())
            {
                break;
            }

            String mismatchDelStr = parseNonIntPrefix(md, i);
            i += mismatchDelStr.length();
            if(mismatchDelStr.charAt(0) == '^')
            {
                for(int j = 1; j < mismatchDelStr.length(); j++)
                {
                    ops.add(new MDOperator(D, mismatchDelStr.charAt(j)));
                }
            }
            else
            {
                ops.add(new MDOperator(M, mismatchDelStr.charAt(0)));
            }
        }

        return ops;
    }

    private static class MismatchCounts
    {
        public final int LowQualMismatchCounts;
        public final int MidQualMismatchCounts;
        public final int HighQualMismatchCounts;

        public MismatchCounts(int lowQualMismatchCounts, int midQualMismatchCounts, int highQualMismatchCounts)
        {
            LowQualMismatchCounts = lowQualMismatchCounts;
            MidQualMismatchCounts = midQualMismatchCounts;
            HighQualMismatchCounts = highQualMismatchCounts;
        }
    }

    private static MismatchCounts getNonMissingMismatches(final SAMRecord read)
    {
        if(read == null || read.getReadUnmappedFlag())
        {
            return new MismatchCounts(0, 0, 0);
        }

        List<MDOperator> mdOperators = parseMDTag(read);
        List<CigarOperator> cigarOps = Lists.newArrayList();
        for(CigarElement el : read.getCigar().getCigarElements())
        {
            for(int i = 0; i < el.getLength(); i++)
            {
                cigarOps.add(el.getOperator());
            }
        }

        int mdIdx = 0;
        int baseIdx = 0;
        String readStr = read.getReadString();
        byte[] baseQuals = read.getBaseQualities();
        int lowQualCount = 0;
        int midQualCount = 0;
        int highQualCount = 0;
        for(CigarOperator op : cigarOps)
        {
            MDOperator mdOp = mdIdx == mdOperators.size() ? null : mdOperators.get(mdIdx);

            if(op == D)
            {
                mdIdx++;
            }
            else if(op == M)
            {
                char base = readStr.charAt(baseIdx);
                int qual = baseQuals[baseIdx];
                if(mdOp.RefBase != null && base != 'N')
                {
                    if(qual <= LOW_QUAL_CUTOFF)
                    {
                        lowQualCount++;
                    }
                    else if(qual <= MID_QUAL_CUTOFF)
                    {
                        midQualCount++;
                    }
                    else
                    {
                        highQualCount++;
                    }
                }

                mdIdx++;
            }

            if(op.consumesReadBases())
            {
                baseIdx++;
            }
        }

        return new MismatchCounts(lowQualCount, midQualCount, highQualCount);
    }

    private void processReadPair(
            final BufferedWriter writer, @Nullable final SAMRecord refRead, @Nullable final SAMRecord newRead) throws IOException
    {
        String readName;
        if(refRead != null)
        {
            readName = refRead.getReadName();
        }
        else
        {
            readName = newRead.getReadName();
        }

        AlignmentInfo refAlignmentInfo = getAlignmentInfo(refRead);
        AlignmentInfo newAlignmentInfo = getAlignmentInfo(newRead);

        SuppInfo refSuppInfo = getSuppInfo(refRead);
        SuppInfo newSuppInfo = getSuppInfo(newRead);

        MismatchCounts refNonMissingMismatches = getNonMissingMismatches(refRead);
        MismatchCounts newNonMissingMismatches = getNonMissingMismatches(newRead);

        String refMDAttribute = refRead == null ? null : refRead.getStringAttribute("MD");
        String newMDAttribute = newRead == null ? null : newRead.getStringAttribute("MD");

        // N bases
        int refNCount = 0;
        if(refRead != null)
        {
            String refReadStr = refRead.getReadString();
            for(int i = 0; i < refReadStr.length(); i++)
            {
                if(refReadStr.charAt(i) == (char) MISSING_BASE)
                {
                    refNCount++;
                }
            }
        }

        int newNCount = 0;
        if(newRead != null)
        {
            String newReadStr = newRead.getReadString();
            for(int i = 0; i < newReadStr.length(); i++)
            {
                if(newReadStr.charAt(i) == (char) MISSING_BASE)
                {
                    newNCount++;
                }
            }
        }

        // write output
        StringJoiner compareLine = new StringJoiner(COMPARE_DELIMITER);
        compareLine.add(readName);
        compareLine.add(refRead == null ? "-" : refRead.getReadString());
        compareLine.add(refRead == null ? "-" : sanatizeQualString(refRead.getBaseQualityString()));

        compareLine.add(newRead == null ? "-" : newRead.getReadString());
        compareLine.add(newRead == null ? "-" : sanatizeQualString(newRead.getBaseQualityString()));

        compareLine.add(refRead == null ? "0" : String.valueOf(refRead.getReadString().length()));
        compareLine.add(refRead == null || refRead.getReadUnmappedFlag() ? "-" : (refRead.getReadNegativeStrandFlag() ? "neg" : "pos"));
        compareLine.add(refRead == null ? "-" : refRead.getCigarString());
        compareLine.add(refRead == null ? "-" : String.valueOf(refRead.getMappingQuality()));

        compareLine.add(refAlignmentInfo == null ? "-" : refAlignmentInfo.Chromosome);
        compareLine.add(refAlignmentInfo == null ? "-" : String.valueOf(refAlignmentInfo.Pos));
        compareLine.add(refAlignmentInfo == null ? "-" : String.valueOf(refAlignmentInfo.AlignmentLength));
        compareLine.add(refRead == null ? "-" : String.valueOf(refNCount));
        compareLine.add(String.valueOf(refNonMissingMismatches.LowQualMismatchCounts));
        compareLine.add(String.valueOf(refNonMissingMismatches.MidQualMismatchCounts));
        compareLine.add(String.valueOf(refNonMissingMismatches.HighQualMismatchCounts));
        compareLine.add(refAlignmentInfo == null ? "-" : String.valueOf(refAlignmentInfo.MatchCount));
        compareLine.add(refAlignmentInfo == null ? "-" : String.valueOf(refAlignmentInfo.LeftSoftClipLength));
        compareLine.add(refAlignmentInfo == null ? "-" : String.valueOf(refAlignmentInfo.RightSoftClipLength));
        compareLine.add(refAlignmentInfo == null
                ? "-"
                : String.valueOf(refAlignmentInfo.LeftSoftClipLength + refAlignmentInfo.RightSoftClipLength));
        compareLine.add(refAlignmentInfo == null ? "-" : String.valueOf(refAlignmentInfo.InsertCount + refAlignmentInfo.DeleteCount));

        compareLine.add(refSuppInfo == null ? "0" : String.valueOf(refSuppInfo.ShortInvSupps.size()));
        compareLine.add(refSuppInfo == null ? "0" : String.valueOf(refSuppInfo.OtherSupps.size()));
        compareLine.add(refSuppInfo == null ? "-" : String.valueOf(refSuppInfo.ShortInvSupps));
        compareLine.add(refSuppInfo == null ? "-" : String.valueOf(refSuppInfo.OtherSupps));
        compareLine.add(refMDAttribute == null ? "-" : refMDAttribute);

        compareLine.add(newRead == null ? "0" : String.valueOf(newRead.getReadString().length()));
        compareLine.add(newRead == null || newRead.getReadUnmappedFlag() ? "-" : (newRead.getReadNegativeStrandFlag() ? "neg" : "pos"));
        compareLine.add(newRead == null ? "-" : newRead.getCigarString());
        compareLine.add(newRead == null ? "-" : String.valueOf(newRead.getMappingQuality()));

        compareLine.add(newAlignmentInfo == null ? "-" : newAlignmentInfo.Chromosome);
        compareLine.add(newAlignmentInfo == null ? "-" : String.valueOf(newAlignmentInfo.Pos));
        compareLine.add(newAlignmentInfo == null ? "-" : String.valueOf(newAlignmentInfo.AlignmentLength));
        compareLine.add(newRead == null ? "-" : String.valueOf(newNCount));
        compareLine.add(String.valueOf(newNonMissingMismatches.LowQualMismatchCounts));
        compareLine.add(String.valueOf(newNonMissingMismatches.MidQualMismatchCounts));
        compareLine.add(String.valueOf(newNonMissingMismatches.HighQualMismatchCounts));
        compareLine.add(newAlignmentInfo == null ? "-" : String.valueOf(newAlignmentInfo.MatchCount));
        compareLine.add(newAlignmentInfo == null ? "-" : String.valueOf(newAlignmentInfo.LeftSoftClipLength));
        compareLine.add(newAlignmentInfo == null ? "-" : String.valueOf(newAlignmentInfo.RightSoftClipLength));
        compareLine.add(newAlignmentInfo == null
                ? "-"
                : String.valueOf(newAlignmentInfo.LeftSoftClipLength + newAlignmentInfo.RightSoftClipLength));
        compareLine.add(newAlignmentInfo == null ? "-" : String.valueOf(newAlignmentInfo.InsertCount + newAlignmentInfo.DeleteCount));

        compareLine.add(newSuppInfo == null ? "0" : String.valueOf(newSuppInfo.ShortInvSupps.size()));
        compareLine.add(newSuppInfo == null ? "0" : String.valueOf(newSuppInfo.OtherSupps.size()));
        compareLine.add(newSuppInfo == null ? "-" : String.valueOf(newSuppInfo.ShortInvSupps));
        compareLine.add(newSuppInfo == null ? "-" : String.valueOf(newSuppInfo.OtherSupps));
        compareLine.add(newMDAttribute == null ? "-" : newMDAttribute);

        writer.write(compareLine.toString());
        writer.newLine();
    }

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder();
        BiomodalCompareConfig.registerConfig(configBuilder);
        configBuilder.checkAndParseCommandLine(args);

        BiomodalCompareConfig config = new BiomodalCompareConfig(configBuilder);
        BiomodalCompare biomodalCompare = new BiomodalCompare(config);
        biomodalCompare.run();
    }
}
