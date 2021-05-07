package com.hartwig.hmftools.lilac.read;

import static com.hartwig.hmftools.lilac.LilacConfig.LL_LOGGER;
import static com.hartwig.hmftools.lilac.fragment.NucleotideFragment.reduceById;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.region.CodingRegions;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.fragment.NucleotideFragment;
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
    private final List<HmfTranscriptRegion> mTranscripts;
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

    // keep track of this common INDEL allele
    private static final Indel STOP_LOSS_ON_C = new Indel("6", 31237115, "CN", "C");

    public SAMRecordReader(
            final String bamFile, final String refGenome, final List<HmfTranscriptRegion> transcripts, final NucleotideFragmentFactory factory)
    {
        mBamFile = bamFile;

        mSamReaderFactory = SamReaderFactory.makeDefault().referenceSequence(new File(refGenome));

        mTranscripts = transcripts;
        mCodingRegions = transcripts.stream()
                .map(x -> new BaseRegion(
                        x.chromosome(), (int)x.codingStart() - MAX_DISTANCE, (int)x.codingEnd() + MAX_DISTANCE))
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

    public final List<NucleotideFragment> readFromBam()
    {
        final List<NucleotideFragment> fragments = Lists.newArrayList();
        mTranscripts.forEach(x -> fragments.addAll(readFromBam(x)));
        return fragments;
    }

    private List<NucleotideFragment> readFromBam(final HmfTranscriptRegion transcript)
    {
        LL_LOGGER.info("  querying {} coding region({}: {} -> {})",
                transcript.gene(), transcript.chromosome(), transcript.codingStart(), transcript.codingEnd());

        boolean reverseStrand = transcript.strand() == Strand.REVERSE;
        final List<NamedBed> codingRegions = codingRegions(transcript, reverseStrand);

        final List<NucleotideFragment> readFragments = Lists.newArrayList();

        for (NamedBed codingRegion : codingRegions)
        {
            readFragments.addAll(realign(codingRegion, reverseStrand, mBamFile));
        }

        return reduceById(readFragments);
    }

    private final List<NamedBed> codingRegions(final HmfTranscriptRegion transcript, boolean reverse)
    {
        final List<NamedBed> regions = CodingRegions.codingRegions(transcript);

        if(!reverse)
            return regions;

        final List<NamedBed> regionsReversed = Lists.newArrayList();

        for(int i = regions.size() - 1; i >= regions.size(); ++i)
            regionsReversed.add(regions.get(i));

        return regionsReversed;
    }

    public final List<NucleotideFragment> readFromBam(final VariantContextDecorator variant)
    {
        return Lists.newArrayList();

        // TODO

        /*
                val variantPosition = GenomePositions.create(variant.chromosome(), variant.position())

        for (transcript in transcripts) {
            val reverseStrand = transcript.strand() == Strand.REVERSE
            val codingRegions = if (reverseStrand) codingRegions(transcript).reversed() else codingRegions(transcript)
            var hlaCodingRegionOffset = 0
            for (codingRegion in codingRegions) {
                if (codingRegion.contains(variantPosition) || abs(codingRegion.start() - variant.position()) <= 5 || abs(codingRegion.end() - variant.position()) <= 5) {
                    val codingRecords = query(reverseStrand, variantPosition, codingRegion, mBamFile)
                            .filter { recordContainsVariant(variant, it) }
                            .distinct()

                    val nucleotideFragments = codingRecords
                            .mapNotNull { factory.createAlignmentFragments(it, codingRegion) }

                    val mateFragments = queryMateFragments(transcript, codingRecords)

                    return (nucleotideFragments + mateFragments)
                            .groupBy { it.id }
                            .map { it.value.reduce { x, y -> NucleotideFragment.merge(x, y) } }

                }
                hlaCodingRegionOffset += codingRegion.bases().toInt()
            }
        }

         */
    }

    private final boolean recordContainsVariant(final VariantContextDecorator variant, final SAMCodingRecord record)
    {
        return false;

        /*
        if(variant.alt().length() != variant.ref().length())
        {
            boolean bl;
            block7:
            {
                String string = variant.chromosome();
                Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "variant.chromosome()");
                int n = (int) variant.position();
                String string2 = variant.ref();
                Intrinsics.checkExpressionValueIsNotNull((Object) string2, (String) "variant.ref()");
                String string3 = variant.alt();
                Intrinsics.checkExpressionValueIsNotNull((Object) string3, (String) "variant.alt()");
                Indel expectedIndel = new Indel(string, n, string2, string3);
                Iterable $receiver$iv = record.getIndels();
                if($receiver$iv instanceof Collection && ((Collection) $receiver$iv).isEmpty())
                {
                    bl = false;
                }
                else
                {
                    for(Object element$iv : $receiver$iv)
                    {
                        Indel it = (Indel) element$iv;
                        boolean bl2 = false;
                        if(!it.match(expectedIndel))
                        {
                            continue;
                        }
                        bl = true;
                        break block7;
                    }
                    bl = false;
                }
            }
            return bl;
        }
        int expectedIndel = 0;
        String string = variant.alt();
        Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "variant.alt()");
        int n = ((CharSequence) string).length();
        while(expectedIndel < n)
        {
            void i;
            int position = (int) variant.position() + i;
            char expectedBase = variant.alt().charAt((int) i);
            int readIndex = record.getRecord().getReadPositionAtReferencePosition(position) - 1;
            if(readIndex < 0)
            {
                return false;
            }
            if((char) record.getRecord().getReadBases()[readIndex] != expectedBase)
            {
                return false;
            }
            ++i;
        }
        return true;

         */
    }

    private final List<NucleotideFragment> queryMateFragments(
            final HmfTranscriptRegion transcript, final List<SAMCodingRecord> codingRecords)
    {
        // TODO for variant calling
        return Lists.newArrayList();

        /*
        val slicer = SAMSlicer(MIN_MAPPING_QUALITY)
        val samRecords = codingRecords.map { it.record }.distinct()
        val mates = samReaderFactory().open(File(mBamFile)).use { reader -> slicer.queryMates(reader, samRecords) }
        val result = mutableListOf<NucleotideFragment>()

        val reverseStrand = transcript.strand() == Strand.REVERSE
        val codingRegions = if (reverseStrand) codingRegions(transcript).reversed() else codingRegions(transcript)

        for (codingRegion in codingRegions) {
            mates
                    .filter { it.alignmentStart <= codingRegion.end() && it.alignmentEnd >= codingRegion.start() }
                    .map { SAMCodingRecord.create(reverseStrand, codingRegion, it) }
                    .mapNotNull { factory.createAlignmentFragments(it, codingRegion) }
                    .forEach { result.add(it) }
        }

        return result

         */
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

    private List<NucleotideFragment> realign(final NamedBed codingRegion, boolean reverseStrand, String bamFileName)
    {
        List<NucleotideFragment> fragments = Lists.newArrayList();

        final List<SAMCodingRecord> codingRecords = query(reverseStrand, codingRegion, bamFileName);

        for(SAMCodingRecord codingRecord : codingRecords)
        {
            if(codingRecord.getIndels().contains(STOP_LOSS_ON_C))
            {
                incrementIndelCounter(mStopLossOnC, STOP_LOSS_ON_C);
            }

            NucleotideFragment fragment = mFragmentFactory.createFragment(codingRecord, codingRegion);

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
