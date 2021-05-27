package com.hartwig.hmftools.lilac.read;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.FusionCommon.NEG_STRAND;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_CHR;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_GENES;
import static com.hartwig.hmftools.lilac.LilacConstants.STOP_LOSS_ON_C;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.LociPosition;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.fragment.FragmentUtils;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragmentFactory;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class SAMRecordReader
{
    private final Map<String,TranscriptData> mTranscripts;
    private final List<BaseRegion> mCodingRegions;

    private final String mBamFile;
    private final SamReaderFactory mSamReaderFactory;
    private final NucleotideFragmentFactory mFragmentFactory;

    private final Map<Indel,Integer> mStopLossOnC;
    private final Map<Indel,Integer> mUnmatchedIndels;
    private final Map<Indel,Integer> mUnmatchedPONIndels;
    private int mFilteredRecordCount;

    public static final int MIN_MAPPING_QUALITY = 1;
    public static final int MAX_DISTANCE = 1000;

    private static final Set<Indel> INDEL_PON = Sets.newHashSet();

    public SAMRecordReader(
            final String bamFile, final String refGenome, final Map<String,TranscriptData> transcripts, final NucleotideFragmentFactory factory)
    {
        mBamFile = bamFile;

        mSamReaderFactory = SamReaderFactory.makeDefault().referenceSequence(new File(refGenome));

        mTranscripts = transcripts;

        mCodingRegions = HLA_GENES.stream()
                .map(x -> mTranscripts.get(x))
                .map(x -> new BaseRegion(HLA_CHR, x.CodingStart - MAX_DISTANCE, x.CodingEnd + MAX_DISTANCE))
                .collect(Collectors.toList());

        mFragmentFactory = factory;
        mFilteredRecordCount = 0;

        mStopLossOnC = Maps.newHashMap();
        mUnmatchedIndels = Maps.newHashMap();
        mUnmatchedPONIndels = Maps.newHashMap();

        // load indel PON
        final List<String> ponLines = new BufferedReader(new InputStreamReader(
                RefGenomeCoordinates.class.getResourceAsStream("/pon/indels.txt")))
                .lines().collect(Collectors.toList());
        ponLines.stream().map(x -> Indel.fromString(x)).forEach(x -> INDEL_PON.add(x));
    }

    public final int alignmentFiltered()
    {
        return mFilteredRecordCount;
    }

    public final int stopLossOnCIndels()
    {
        return mStopLossOnC.containsKey(STOP_LOSS_ON_C) ? mStopLossOnC.get(STOP_LOSS_ON_C) : 0;
    }

    public final Map<Indel,Integer> unmatchedIndels(int minCount)
    {
        Map<Indel,Integer> filteredMap = Maps.newHashMap();
        mUnmatchedIndels.entrySet().stream().filter(x -> x.getValue()>= minCount).forEach(x -> filteredMap.put(x.getKey(), x.getValue()));
        return filteredMap;
    }

    public final Map<Indel, Integer> unmatchedPonIndels(int minCount)
    {
        Map<Indel,Integer> filteredMap = Maps.newHashMap();
        mUnmatchedPONIndels.entrySet().stream().filter(x -> x.getValue()>= minCount).forEach(x -> filteredMap.put(x.getKey(), x.getValue()));
        return filteredMap;
    }

    public final List<Fragment> readFromBam()
    {
        final List<Fragment> fragments = Lists.newArrayList();

        for(String geneName : HLA_GENES)
        {
            TranscriptData transcript = mTranscripts.get(geneName);
            fragments.addAll(readFromBam(geneName, transcript));
        }

        return fragments;
    }

    private List<Fragment> readFromBam(final String geneName, final TranscriptData transcript)
    {
        LL_LOGGER.info("  querying {} coding region({}: {} -> {})",
                geneName, HLA_CHR, transcript.CodingStart, transcript.CodingEnd);

        boolean reverseStrand = transcript.Strand == NEG_STRAND;
        final List<NamedBed> codingRegions = codingRegions(geneName, transcript);

        final List<Fragment> readFragments = Lists.newArrayList();

        for (NamedBed codingRegion : codingRegions)
        {
            readFragments.addAll(realign(codingRegion, reverseStrand, mBamFile));
        }

        return FragmentUtils.mergeFragmentsById(readFragments);
    }

    private final List<NamedBed> codingRegions(final String geneName, final TranscriptData transcript)
    {
        final List<NamedBed> regions = LociPosition.codingRegions(geneName, HLA_CHR, transcript);

        if(transcript.Strand == POS_STRAND)
            return regions;

        final List<NamedBed> regionsReversed = Lists.newArrayList();

        for(int i = regions.size() - 1; i >= 0; --i)
            regionsReversed.add(regions.get(i));

        return regionsReversed;
    }

    public final List<Fragment> readFromBam(final VariantContextDecorator variant)
    {
        final GenomePosition variantPosition = GenomePositions.create(variant.chromosome(), variant.position());

        for(String geneName : HLA_GENES)
        {
            TranscriptData transcript = mTranscripts.get(geneName);
            boolean reverseStrand = transcript.Strand == NEG_STRAND;
            final List<NamedBed> codingRegions = codingRegions(geneName, transcript);

            for(NamedBed codingRegion : codingRegions)
            {
                if (codingRegion.contains(variantPosition)
                || abs(codingRegion.start() - variant.position()) <= 5 || abs(codingRegion.end() - variant.position()) <= 5)
                {
                    List<SAMCodingRecord> regionCodingRecords = query(reverseStrand, variantPosition, codingRegion, mBamFile);
                    List<SAMCodingRecord> codingRecords = Lists.newArrayList();

                    for(SAMCodingRecord record : regionCodingRecords)
                    {
                        if(!recordContainsVariant(variant, record))
                            continue;

                        if(codingRecords.stream().anyMatch(x -> x.getSamRecord().hashCode() == record.getSamRecord().hashCode()))
                            continue;

                        codingRecords.add(record);
                    }

                    final List<Fragment> readFragments = codingRecords.stream()
                            .map(x -> mFragmentFactory.createAlignmentFragments(x, codingRegion))
                            .filter(x -> x != null)
                            .collect(Collectors.toList());

                    List<Fragment> mateFragments = queryMateFragments(geneName, transcript, codingRecords);
                    readFragments.addAll(mateFragments);

                    List<Fragment> readGroupFragments = FragmentUtils.mergeFragmentsById(readFragments);
                    return readGroupFragments;
                }
            }
        }

        return Lists.newArrayList();
    }

    private final boolean recordContainsVariant(final VariantContextDecorator variant, final SAMCodingRecord record)
    {
        if (variant.alt().length() != variant.ref().length())
        {
            Indel expectedIndel = new Indel(variant.chromosome(), (int)variant.position(), variant.ref(), variant.alt());
            return record.getIndels().stream().anyMatch(x -> x.match(expectedIndel));
        }

        for (int i = 0; i < variant.alt().length(); ++i)
        {
            int position = (int)variant.position() + i;
            char expectedBase = variant.alt().charAt(i);

            int readIndex = record.getSamRecord().getReadPositionAtReferencePosition(position) - 1;

            if (readIndex < 0)
                return false;

            if (record.getSamRecord().getReadString().charAt(readIndex) != expectedBase)
                return false;
        }

        return true;
    }

    private final List<Fragment> queryMateFragments(
            final String geneName, final TranscriptData transcript, final List<SAMCodingRecord> codingRecords)
    {
        SAMSlicer slicer = new SAMSlicer(MIN_MAPPING_QUALITY);

        SamReader samReader = mSamReaderFactory.open(new File(mBamFile));

        List<SAMRecord> records = codingRecords.stream().map(x -> x.getSamRecord()).collect(Collectors.toList());
        List<SAMRecord> mateRecords = slicer.queryMates(samReader, records);

        List<Fragment> fragments = Lists.newArrayList();

        boolean reverseStrand = transcript.Strand == NEG_STRAND;
        final List<NamedBed> codingRegions = codingRegions(geneName, transcript);

        for (NamedBed codingRegion : codingRegions)
        {
            for(SAMRecord record : mateRecords)
            {
                if(record.getAlignmentStart() > codingRegion.end() || record.getAlignmentEnd() < codingRegion.start())
                    continue;

                SAMCodingRecord codingRecord = SAMCodingRecord.create(
                        reverseStrand, BaseRegion.from(codingRegion), record, true, true);

                Fragment fragment = mFragmentFactory.createAlignmentFragments(codingRecord, codingRegion);

                if(fragment != null)
                    fragments.add(fragment);
            }
        }

        return fragments;
    }

    private List<SAMCodingRecord> query(
            boolean reverseStrand, final GenomePosition variantRegion, final NamedBed nearestCodingRegion, final String bamFileName)
    {
        SAMSlicer slicer = new SAMSlicer(MIN_MAPPING_QUALITY);

        SamReader samReader = mSamReaderFactory.open(new File(bamFileName));

        BaseRegion codingRegion = new BaseRegion(
                nearestCodingRegion.chromosome(), (int)nearestCodingRegion.start(), (int)nearestCodingRegion.end());

        final List<SAMRecord> records = slicer.slice(
                variantRegion.chromosome(), (int)variantRegion.position(), (int)variantRegion.position(), samReader);

        final List<SAMCodingRecord> codingRecords = Lists.newArrayList();

        for(SAMRecord record : records)
        {
            if(bothEndsInRangeOfCodingTranscripts(record))
            {
                codingRecords.add(SAMCodingRecord.create(reverseStrand, codingRegion, record, true, true));
            }
            else
            {
                ++mFilteredRecordCount;
            }
        }

        return codingRecords;
    }

    private List<SAMCodingRecord> query(boolean reverseStrand, final NamedBed bedRegion, final String bamFileName)
    {
        SAMSlicer slicer = new SAMSlicer(MIN_MAPPING_QUALITY);

        SamReader samReader = mSamReaderFactory.open(new File(bamFileName));

        BaseRegion codingRegion = new BaseRegion(bedRegion.chromosome(), (int)bedRegion.start(), (int)bedRegion.end());

        final List<SAMRecord> records = slicer.slice(codingRegion.Chromosome, codingRegion.start(), codingRegion.end(), samReader);

        final List<SAMCodingRecord> codingRecords = Lists.newArrayList();

        for(SAMRecord record : records)
        {
            if(bothEndsInRangeOfCodingTranscripts(record))
            {
                codingRecords.add(SAMCodingRecord.create(reverseStrand, codingRegion, record, true, true));
            }
            else
            {
                ++mFilteredRecordCount;
            }
        }

        return codingRecords;
    }

    private List<Fragment> realign(final NamedBed codingRegion, boolean reverseStrand, String bamFileName)
    {
        List<Fragment> fragments = Lists.newArrayList();

        final List<SAMCodingRecord> codingRecords = query(reverseStrand, codingRegion, bamFileName);

        for(SAMCodingRecord codingRecord : codingRecords)
        {
            if(codingRecord.getIndels().contains(STOP_LOSS_ON_C))
            {
                incrementIndelCounter(mStopLossOnC, STOP_LOSS_ON_C);
            }

            Fragment fragment = mFragmentFactory.createFragment(codingRecord, codingRegion);

            if(fragment != null)
            {
                fragments.add(fragment);
                continue;
            }

            for(Indel indel : codingRecord.getIndels())
            {
                if(INDEL_PON.contains(indel))
                {
                    incrementIndelCounter(mUnmatchedPONIndels, indel);
                    continue;
                }

                incrementIndelCounter(mUnmatchedIndels, indel);
            }
        }

        return fragments;
    }

    private static void incrementIndelCounter(final Map<Indel,Integer> indelMap, final Indel indel)
    {
        Integer count = indelMap.get(indel);
        indelMap.put(indel, count != null ? count + 1 : 1);
    }

    private final boolean bothEndsInRangeOfCodingTranscripts(final SAMRecord record)
    {
        boolean recordInRange = mCodingRegions.stream()
                .anyMatch(x -> x.containsPosition(record.getContig(), record.getAlignmentStart()));

        boolean mateInRange = mCodingRegions.stream()
                .anyMatch(x -> x.containsPosition(record.getContig(), record.getMateAlignmentStart()));

        return recordInRange && mateInRange;
    }
}
