package com.hartwig.hmftools.lilac.read;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_CHR;
import static com.hartwig.hmftools.lilac.LilacConstants.HLA_GENES;
import static com.hartwig.hmftools.lilac.LilacConstants.SPLICE_VARIANT_BUFFER;
import static com.hartwig.hmftools.lilac.ReferenceData.STOP_LOSS_ON_C_INDEL;
import static com.hartwig.hmftools.lilac.ReferenceData.INDEL_PON;
import static com.hartwig.hmftools.lilac.fragment.FragmentUtils.mergeFragments;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.lilac.LilacConfig;
import com.hartwig.hmftools.lilac.fragment.Fragment;
import com.hartwig.hmftools.lilac.fragment.FragmentUtils;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragmentFactory;
import com.hartwig.hmftools.lilac.variant.SomaticVariant;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class BamRecordReader implements BamReader
{
    private final Map<String,GeneCodingRegions> mGeneCodingRegions;

    private final String mBamFile;
    private final SamReader mSamReader;

    private final BamSlicer mBamSlicer;
    private final NucleotideFragmentFactory mFragmentFactory;

    private final Map<Indel,List<Fragment>> mKnownStopLossFragments;
    private final Map<Indel,Integer> mUnmatchedIndels;
    private final Map<Indel,Integer> mUnmatchedPONIndels;
    private final Set<String> mDiscardIndelReadIds;
    private int mFilteredRecordCount;

    public static final int MIN_MAPPING_QUALITY = 0;
    public static final int MAX_DISTANCE = 1000;

    public BamRecordReader(
            final String bamFile, final LilacConfig config, final Map<String,TranscriptData> transcripts, final NucleotideFragmentFactory factory)
    {
        mBamFile = bamFile;

        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault()
                .validationStringency(config.BamStringency)
                .referenceSequence(new File(config.RefGenome));

        mSamReader = samReaderFactory.open(new File(mBamFile));
        mBamSlicer = new BamSlicer(MIN_MAPPING_QUALITY);

        mGeneCodingRegions = Maps.newHashMap();

        for(String geneName : HLA_GENES)
        {
            TranscriptData transcriptData = transcripts.get(geneName);
            mGeneCodingRegions.put(geneName, new GeneCodingRegions(geneName, HLA_CHR, transcriptData));
        }

        mFragmentFactory = factory;
        mFilteredRecordCount = 0;

        mKnownStopLossFragments = Maps.newHashMap();
        mUnmatchedIndels = Maps.newHashMap();
        mUnmatchedPONIndels = Maps.newHashMap();
        mDiscardIndelReadIds = Sets.newHashSet();
    }

    public int filteredReadCount() { return mFilteredRecordCount; }

    public Map<Indel,List<Fragment>> getKnownStopLossFragments() { return mKnownStopLossFragments; }

    public Map<Indel,Integer> unmatchedIndels(int minCount)
    {
        Map<Indel,Integer> filteredMap = Maps.newHashMap();
        mUnmatchedIndels.entrySet().stream()
                .filter(x -> x.getValue()>= minCount)
                .filter(x -> !x.getKey().equals(STOP_LOSS_ON_C_INDEL))
                .forEach(x -> filteredMap.put(x.getKey(), x.getValue()));
        return filteredMap;
    }

    public Map<Indel,Integer> unmatchedPonIndels(int minCount)
    {
        Map<Indel,Integer> filteredMap = Maps.newHashMap();
        mUnmatchedPONIndels.entrySet().stream().filter(x -> x.getValue()>= minCount).forEach(x -> filteredMap.put(x.getKey(), x.getValue()));
        return filteredMap;
    }

    public List<Fragment> findGeneFragments()
    {
        final List<Fragment> fragments = Lists.newArrayList();

        for(String geneName : HLA_GENES)
        {
            GeneCodingRegions geneCodingRegions = mGeneCodingRegions.get(geneName);

            fragments.addAll(findGeneFragments(geneCodingRegions));
        }

        return fragments;
    }

    private List<Fragment> findGeneFragments(final GeneCodingRegions geneCodingRegions)
    {
        LL_LOGGER.trace("querying HLA gene({})", geneCodingRegions.GeneName);

        // slice for the whole coding region rather than per exon
        ChrBaseRegion sliceRegion = new ChrBaseRegion(geneCodingRegions.Chromosome, geneCodingRegions.CodingStart, geneCodingRegions.CodingEnd);
        List<SAMRecord> records = mBamSlicer.slice(mSamReader, sliceRegion);

        List<ReadRecord> reads = Lists.newArrayList();

        for(SAMRecord record : records)
        {
            List<BaseRegion> codingExonOverlaps = geneCodingRegions.CodingRegions.stream()
                    .filter(x -> positionsOverlap(x.start(), x.end(), record.getAlignmentStart(), record.getAlignmentEnd()))
                    .collect(Collectors.toList());

            // read must overlap an exon to even be considered
            if(codingExonOverlaps.isEmpty())
                continue;

            if(!bothEndsInRangeOfCodingTranscripts(record))
            {
                LL_LOGGER.trace("filter read: id({}) coords({}-{}) cigar({}) from gene region({}:{}-{})",
                        record.getReadName(), record.getAlignmentStart(), record.getAlignmentEnd(), record.getCigarString(),
                        geneCodingRegions.GeneName, geneCodingRegions.CodingStart, geneCodingRegions.CodingEnd);

                ++mFilteredRecordCount;
            }
            else
            {
                codingExonOverlaps.forEach(x -> reads.add(ReadRecord.create(x, record, true, true)));
            }
        }

        Map<String,Fragment> readIdFragments = createFragments(geneCodingRegions.GeneName, geneCodingRegions.Strand, reads);

        List<Fragment> readFragments = readIdFragments.values().stream()
                .filter(x -> !mDiscardIndelReadIds.contains(x.id()))
                .collect(Collectors.toList());

        return readFragments;
    }

    private Map<String,Fragment> createFragments(final String geneName, final byte geneStrand, final List<ReadRecord> codingRecords)
    {
        Map<String,Fragment> readIdFragments = Maps.newHashMap();

        for(ReadRecord codingRecord : codingRecords)
        {
            Fragment fragment = mFragmentFactory.createFragment(codingRecord, geneName, geneStrand);

            if(fragment != null)
            {
                Fragment existingFragment = readIdFragments.get(fragment.id());

                if(existingFragment == null)
                {
                    readIdFragments.put(fragment.id(), fragment);
                    existingFragment = fragment;
                }
                else
                {
                    mergeFragments(existingFragment, fragment);
                }

                if(codingRecord.getIndels().contains(STOP_LOSS_ON_C_INDEL))
                    addKnownIndelFragment(existingFragment);

                continue;
            }

            if(codingRecord.containsIndel())
            {
                if(codingRecord.getIndels().contains(STOP_LOSS_ON_C_INDEL))
                {
                    LL_LOGGER.trace("missing known indel fragment: {} {}", codingRecord.Id, codingRecord.readInfo());
                }

                // the other read belonging to this fragment won't be used
                mDiscardIndelReadIds.add(codingRecord.Id);
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

        return readIdFragments;
    }

    private void addKnownIndelFragment(final Fragment fragment)
    {
        List<Fragment> indelFrags = mKnownStopLossFragments.get(STOP_LOSS_ON_C_INDEL);
        if(indelFrags == null)
        {
            indelFrags = Lists.newArrayList();
            mKnownStopLossFragments.put(STOP_LOSS_ON_C_INDEL, indelFrags);
        }
        else if(indelFrags.contains(fragment))
        {
            return;
        }

        indelFrags.add(fragment);
    }

    private static void incrementIndelCounter(final Map<Indel,Integer> indelMap, final Indel indel)
    {
        Integer count = indelMap.get(indel);
        indelMap.put(indel, count != null ? count + 1 : 1);
    }

    private boolean bothEndsInRangeOfCodingTranscripts(final SAMRecord record)
    {
        if(!record.getMateReferenceName().equals(HLA_CHR))
            return false;

        // this check allows records to span across HLA genes, since a read may be mismapped
        boolean readInRange = mGeneCodingRegions.values().stream()
                .anyMatch(x -> x.withinCodingBounds(record.getAlignmentStart(), MAX_DISTANCE));

        boolean mateInRange = mGeneCodingRegions.values().stream()
                .anyMatch(x -> x.withinCodingBounds(record.getMateAlignmentStart(), MAX_DISTANCE));

        return readInRange && mateInRange;
    }

    // methods for checking somatic variant support in HLA regions
    public List<Fragment> findVariantFragments(final SomaticVariant variant)
    {
        // slice the BAM for the variant, get the mates if they are within the same gene's coding region,
        // then create and filter fragments

        for(String geneName : HLA_GENES)
        {
            GeneCodingRegions geneCodingRegions = mGeneCodingRegions.get(geneName);

            for(BaseRegion codingRegion : geneCodingRegions.CodingRegions)
            {
                if(!regionCoversVariant(codingRegion, variant))
                    continue;

                ChrBaseRegion sliceRegion = new ChrBaseRegion(geneCodingRegions.Chromosome, codingRegion.start(), codingRegion.end());
                List<ReadRecord> regionCodingRecords = findVariantRecords(variant, sliceRegion);
                List<ReadRecord> codingRecords = Lists.newArrayList();

                for(ReadRecord record : regionCodingRecords)
                {
                    if(!recordContainsVariant(variant, record))
                        continue;

                    if(codingRecords.stream().anyMatch(x -> x.getSamRecord().hashCode() == record.getSamRecord().hashCode()))
                        continue;

                    codingRecords.add(record);
                }

                final List<Fragment> readFragments = codingRecords.stream()
                        .map(x -> mFragmentFactory.createAlignmentFragments(x, geneName, geneCodingRegions.Strand))
                        .filter(x -> x != null)
                        .collect(Collectors.toList());

                List<Fragment> mateFragments = queryMateFragments(geneCodingRegions, codingRecords);
                readFragments.addAll(mateFragments);

                List<Fragment> readGroupFragments = FragmentUtils.mergeFragmentsById(readFragments);
                return readGroupFragments;
            }
        }

        return Lists.newArrayList();
    }

    private List<Fragment> queryMateFragments(GeneCodingRegions geneCodingRegions, final List<ReadRecord> codingRecords)
    {
        List<SAMRecord> records = codingRecords.stream().map(x -> x.getSamRecord()).collect(Collectors.toList());
        List<SAMRecord> mateRecords = mBamSlicer.queryMates(mSamReader, records);

        List<Fragment> fragments = Lists.newArrayList();

        for(SAMRecord record : mateRecords)
        {
            // take any mate within a coding region of the same gene
            BaseRegion matchedCodingRegion = geneCodingRegions.CodingRegions.stream()
                    .filter(x -> positionsOverlap(x.start(), x.end(), record.getAlignmentStart(), record.getAlignmentEnd()))
                    .findFirst().orElse(null);

            if(matchedCodingRegion == null)
                continue;

            ReadRecord codingRecord = ReadRecord.create(matchedCodingRegion, record, true, true);

            Fragment fragment = mFragmentFactory.createAlignmentFragments(codingRecord, geneCodingRegions.GeneName, geneCodingRegions.Strand);

            if(fragment != null)
                fragments.add(fragment);
        }

        return fragments;
    }

    private List<ReadRecord> findVariantRecords(final SomaticVariant variant, final ChrBaseRegion codingRegion)
    {
        final List<SAMRecord> records = mBamSlicer.slice(mSamReader, new ChrBaseRegion(variant.Chromosome, variant.Position, variant.Position));

        final List<ReadRecord> codingRecords = Lists.newArrayList();

        BaseRegion baseCodingRegion = new BaseRegion(codingRegion.start(), codingRegion.end());

        for(SAMRecord record : records)
        {
            if(bothEndsInRangeOfCodingTranscripts(record))
            {
                codingRecords.add(ReadRecord.create(baseCodingRegion, record, true, true));
            }
            else
            {
                ++mFilteredRecordCount;
            }
        }

        return codingRecords;
    }

    private static boolean regionCoversVariant(final BaseRegion codingRegion, final SomaticVariant variant)
    {
        // allow a margin for splice variants
        if(variant.CanonicalCodingEffect == CodingEffect.SPLICE)
        {
            return abs(variant.Position - codingRegion.start()) <= SPLICE_VARIANT_BUFFER
                    || abs(variant.Position - codingRegion.end()) <= SPLICE_VARIANT_BUFFER;
        }
        else
        {
            return codingRegion.containsPosition(variant.Position);
        }
    }

    public static List<Fragment> filterVariantFragments(final SomaticVariant variant, final List<Fragment> fragments)
    {
        return fragments.stream()
                .filter(x -> x.reads().stream().anyMatch(y -> recordContainsVariant(variant, y)))
                .collect(Collectors.toList());
    }

    private static boolean recordContainsVariant(final SomaticVariant variant, final ReadRecord record)
    {
        if(variant.Alt.length() != variant.Ref.length())
        {
            Indel expectedIndel = new Indel(variant.Chromosome, variant.Position, variant.Ref, variant.Alt);
            return record.getIndels().stream().anyMatch(x -> x.match(expectedIndel));
        }

        for(int i = 0; i < variant.Alt.length(); ++i)
        {
            int position = variant.Position + i;
            char expectedBase = variant.Alt.charAt(i);

            int readIndex = record.getSamRecord().getReadPositionAtReferencePosition(position) - 1;

            if(readIndex < 0)
                return false;

            if(record.getSamRecord().getReadString().charAt(readIndex) != expectedBase)
                return false;
        }

        return true;
    }
}
