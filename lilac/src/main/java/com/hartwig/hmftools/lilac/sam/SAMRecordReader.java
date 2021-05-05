package com.hartwig.hmftools.lilac.sam;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeCoordinates;
import com.hartwig.hmftools.common.genome.region.CodingRegions;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.utils.sv.BaseRegion;
import com.hartwig.hmftools.common.variant.VariantContextDecorator;
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment;
import com.hartwig.hmftools.lilac.nuc.NucleotideFragmentFactory;
import com.hartwig.hmftools.lilac.sam.Indel;
import com.hartwig.hmftools.lilac.sam.SAMCodingRecord;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Consumer;
import java.util.stream.Collectors;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public class SAMRecordReader
{
    private final List<BaseRegion> mCodingRegions;

    private final String mBamFile;
    private final String mRefGenome;
    private final NucleotideFragmentFactory mFragmentFactory;

    private final Map<Indel,Integer> mStopLossOnC;
    private final Map<Indel,Integer> mUnmatchedIndels;
    private final Map<Indel,Integer> mUnmatchedPONIndels;
    private int mFilteredRecordCount;

    public static final int MIN_MAPPING_QUALITY = 1;
    public static final int MAX_DISTANCE = 1000;

    private static final Set<Indel> INDEL_PON = Sets.newHashSet();

    private static final Indel STOP_LOSS_ON_C = new Indel("6", 31237115, "CN", "C");

    public SAMRecordReader(
            final String bamFile, final String refGenome, final List<HmfTranscriptRegion> transcripts, final NucleotideFragmentFactory factory)
    {
        mBamFile = bamFile;
        mRefGenome = refGenome;

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
        Integer n = mStopLossOnC.get(STOP_LOSS_ON_C);
        return n != null ? n : 0;
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
        return Lists.newArrayList();
        // mCodingRegions.forEach(x -> readFromBam());
    }

    private final List<NucleotideFragment> readFromBam(HmfTranscriptRegion transcript, String bamFile)
    {
        return Lists.newArrayList();

        /*
        Object list$iv$iv;
        Map $receiver$iv$iv;
        Map $receiver$iv;
        Iterable codingRegion2;
        logger.info("    querying " + transcript.gene() + " (" + transcript.chromosome() + ':' + transcript.codingStart() + '-' + transcript
                .codingEnd() + ')');
        boolean reverseStrand = transcript.strand() == Strand.REVERSE;
        List codingRegions =
                reverseStrand ? CollectionsKt.reversed((Iterable) codingRegions(transcript)) : codingRegions(transcript);
        List realignedRegions = new ArrayList();
        for(Iterable codingRegion2 : codingRegions)
        {
            realignedRegions.addAll((Collection) realign((NamedBed) codingRegion2, reverseStrand, bamFile));
        }
        codingRegion2 = realignedRegions;
        Iterator iterator = $receiver$iv;
        Object destination$iv$iv = new LinkedHashMap();
        Object object = $receiver$iv$iv.iterator();
        while(object.hasNext())
        {
            Object object2;
            Object element$iv$iv = object.next();
            NucleotideFragment it = (NucleotideFragment) element$iv$iv;
            boolean bl = false;
            Object $receiver$iv$iv$iv = destination$iv$iv;
            String key$iv$iv = it.getId();
            Object value$iv$iv$iv = $receiver$iv$iv$iv.get(key$iv$iv);
            if(value$iv$iv$iv == null)
            {
                ArrayList answer$iv$iv$iv = new ArrayList();
                $receiver$iv$iv$iv.put(key$iv$iv, answer$iv$iv$iv);
                object2 = answer$iv$iv$iv;
            }
            else
            {
                object2 = value$iv$iv$iv;
            }
            list$iv$iv = (List) object2;
            list$iv$iv.add(element$iv$iv);
        }
        $receiver$iv = destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList($receiver$iv.size());
        object = $receiver$iv$iv;
        Iterator iterator2 = object.entrySet().iterator();
        while(iterator2.hasNext())
        {
            void it;
            Map.Entry item$iv$iv;
            Map.Entry bl = item$iv$iv = iterator2.next();
            Object object3 = destination$iv$iv;
            boolean bl2 = false;
            Iterable $receiver$iv2 = (Iterable) it.getValue();
            Iterator iterator$iv = $receiver$iv2.iterator();
            if(!iterator$iv.hasNext())
            {
                throw (Throwable) new UnsupportedOperationException("Empty collection can't be reduced.");
            }
            Object accumulator$iv = iterator$iv.next();
            while(iterator$iv.hasNext())
            {
                void y;
                list$iv$iv = (NucleotideFragment) iterator$iv.next();
                NucleotideFragment x = (NucleotideFragment) accumulator$iv;
                boolean bl3 = false;
                accumulator$iv = NucleotideFragment.Companion.merge(x, (NucleotideFragment) y);
            }
            NucleotideFragment nucleotideFragment = (NucleotideFragment) accumulator$iv;
            object3.add(nucleotideFragment);
        }
        return (List) destination$iv$iv;

         */
    }

    public final List<NucleotideFragment> readFromBam(final VariantContextDecorator variant)
    {
        return Lists.newArrayList();

        /*
                val variantPosition = GenomePositions.create(variant.chromosome(), variant.position())

        for (transcript in transcripts) {
            val reverseStrand = transcript.strand() == Strand.REVERSE
            val codingRegions = if (reverseStrand) codingRegions(transcript).reversed() else codingRegions(transcript)
            var hlaCodingRegionOffset = 0
            for (codingRegion in codingRegions) {
                if (codingRegion.contains(variantPosition) || abs(codingRegion.start() - variant.position()) <= 5 || abs(codingRegion.end() - variant.position()) <= 5) {
                    val codingRecords = query(reverseStrand, variantPosition, codingRegion, bamFile)
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

        /*
        Intrinsics.checkParameterIsNotNull((Object) variant, (String) "variant");
        GenomePosition variantPosition = GenomePositions.create((String) variant.chromosome(), (long) variant.position());
        for(HmfTranscriptRegion transcript : transcripts)
        {
            boolean reverseStrand = transcript.strand() == Strand.REVERSE;
            List codingRegions =
                    reverseStrand ? CollectionsKt.reversed((Iterable) codingRegions(transcript)) : codingRegions(transcript);
            int hlaCodingRegionOffset = 0;
            for(NamedBed codingRegion : codingRegions)
            {
                long l;
                if(codingRegion.contains(variantPosition) || Math.abs(l = codingRegion.start() - variant.position()) <= (long) 5
                        || Math.abs(l = codingRegion.end() - variant.position()) <= (long) 5)
                {
                    Object list$iv$iv;
                    Map $receiver$iv$iv;
                    Map $receiver$iv;
                    Iterable $receiver$iv$iv2;
                    void $receiver$iv2;
                    Iterable $receiver$iv$iv3;
                    void $receiver$iv3;
                    GenomePosition genomePosition = variantPosition;
                    Intrinsics.checkExpressionValueIsNotNull((Object) genomePosition, (String) "variantPosition");
                    Iterable iterable = query(reverseStrand, genomePosition, codingRegion, mBamFile);
                    void var12_13 = $receiver$iv3;
                    Collection destination$iv$iv = new ArrayList();
                    for(Object element$iv$iv : $receiver$iv$iv3)
                    {
                        SAMCodingRecord it = (SAMCodingRecord) element$iv$iv;
                        boolean bl = false;
                        if(!recordContainsVariant(variant, it))
                        {
                            continue;
                        }
                        destination$iv$iv.add(element$iv$iv);
                    }
                    List codingRecords = CollectionsKt.distinct((Iterable) ((List) destination$iv$iv));
                    $receiver$iv$iv3 = codingRecords;
                    destination$iv$iv = $receiver$iv2;
                    Collection destination$iv$iv2 = new ArrayList();
                    void $receiver$iv$iv$iv = $receiver$iv$iv2;
                    Object object = $receiver$iv$iv$iv.iterator();
                    while(object.hasNext())
                    {
                        NucleotideFragment nucleotideFragment;
                        Object element$iv$iv$iv;
                        Object element$iv$iv = element$iv$iv$iv = object.next();
                        SAMCodingRecord it = (SAMCodingRecord) element$iv$iv;
                        boolean bl = false;
                        if(mFragmentFactory.createAlignmentFragments(it, codingRegion) == null)
                        {
                            continue;
                        }
                        NucleotideFragment it$iv$iv = nucleotideFragment;
                        destination$iv$iv2.add(it$iv$iv);
                    }
                    List nucleotideFragments = (List) destination$iv$iv2;
                    List<NucleotideFragment> mateFragments = queryMateFragments(transcript, codingRecords);
                    $receiver$iv$iv2 = CollectionsKt.plus((Collection) nucleotideFragments, (Iterable) mateFragments);
                    destination$iv$iv2 = $receiver$iv;
                    Object destination$iv$iv3 = new LinkedHashMap();
                    object = $receiver$iv$iv.iterator();
                    while(object.hasNext())
                    {
                        Object object2;
                        Object element$iv$iv = object.next();
                        NucleotideFragment it = (NucleotideFragment) element$iv$iv;
                        boolean bl = false;
                        Object $receiver$iv$iv$iv2 = destination$iv$iv3;
                        String key$iv$iv = it.getId();
                        Object value$iv$iv$iv = $receiver$iv$iv$iv2.get(key$iv$iv);
                        if(value$iv$iv$iv == null)
                        {
                            ArrayList answer$iv$iv$iv = new ArrayList();
                            $receiver$iv$iv$iv2.put(key$iv$iv, answer$iv$iv$iv);
                            object2 = answer$iv$iv$iv;
                        }
                        else
                        {
                            object2 = value$iv$iv$iv;
                        }
                        list$iv$iv = (List) object2;
                        list$iv$iv.add(element$iv$iv);
                    }
                    $receiver$iv = destination$iv$iv3;
                    $receiver$iv$iv = $receiver$iv;
                    destination$iv$iv3 = new ArrayList($receiver$iv.size());
                    object = $receiver$iv$iv;
                    Iterator iterator = object.entrySet().iterator();
                    while(iterator.hasNext())
                    {
                        void it;
                        Map.Entry item$iv$iv;
                        Map.Entry bl = item$iv$iv = iterator.next();
                        Object object3 = destination$iv$iv3;
                        boolean bl2 = false;
                        Iterable $receiver$iv4 = (Iterable) it.getValue();
                        Iterator iterator$iv = $receiver$iv4.iterator();
                        if(!iterator$iv.hasNext())
                        {
                            throw (Throwable) new UnsupportedOperationException("Empty collection can't be reduced.");
                        }
                        Object accumulator$iv = iterator$iv.next();
                        while(iterator$iv.hasNext())
                        {
                            void y;
                            list$iv$iv = (NucleotideFragment) iterator$iv.next();
                            NucleotideFragment x = (NucleotideFragment) accumulator$iv;
                            boolean bl3 = false;
                            accumulator$iv = NucleotideFragment.Companion.merge(x, (NucleotideFragment) y);
                        }
                        NucleotideFragment nucleotideFragment = (NucleotideFragment) accumulator$iv;
                        object3.add(nucleotideFragment);
                    }
                    return (List) destination$iv$iv3;
                }
                hlaCodingRegionOffset += (int) codingRegion.bases();
            }
        }
        return CollectionsKt.emptyList();

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


    private final List<NamedBed> codingRegions(HmfTranscriptRegion transcript)
    {
        return CodingRegions.codingRegions((HmfTranscriptRegion) transcript);
    }

    private final SamReaderFactory samReaderFactory()
    {
        SamReaderFactory samReaderFactory;
        SamReaderFactory samReaderFactory2 = SamReaderFactory.makeDefault();
        CharSequence charSequence = mRefGenome;
        if(charSequence.length() > 0)
        {
            SamReaderFactory samReaderFactory3 = samReaderFactory2.referenceSequence(new File(mRefGenome));
            samReaderFactory = samReaderFactory3;
        }
        else
        {
            SamReaderFactory samReaderFactory4 = samReaderFactory2;
            samReaderFactory = samReaderFactory4;
        }
        return samReaderFactory;
    }

    private final List<NucleotideFragment> queryMateFragments(HmfTranscriptRegion transcript, List<SAMCodingRecord> codingRecords)
    {
        return Lists.newArrayList();

        /*
        Object object;
        Collection collection;
        Closeable $receiver$iv$iv;
        Iterable $receiver$iv;
        SAMSlicer slicer = new SAMSlicer(1);
        Iterable iterable = $receiver$iv = (Iterable) codingRecords;
        Object destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        Iterator iterator = $receiver$iv$iv.iterator();
        while(iterator.hasNext())
        {
            void it;
            Object item$iv$iv = iterator.next();
            SAMCodingRecord sAMCodingRecord = (SAMCodingRecord) item$iv$iv;
            collection = destination$iv$iv;
            boolean bl = false;
            object = it.getRecord();
            collection.add(object);
        }
        List samRecords = CollectionsKt.distinct((Iterable) ((List) destination$iv$iv));
        $receiver$iv$iv = (Closeable) samReaderFactory().open(new File(mBamFile));
        destination$iv$iv = null;
        try
        {
            Object reader = (SamReader) $receiver$iv$iv;
            boolean bl = false;
            SamReader samReader = reader;
            Intrinsics.checkExpressionValueIsNotNull((Object) samReader, (String) "reader");
            reader = slicer.queryMates(samReader, samRecords);
        } catch (Throwable reader)
        {
            destination$iv$iv = reader;
            throw reader;
        } finally
        {
            CloseableKt.closeFinally((Closeable) $receiver$iv$iv, (Throwable) destination$iv$iv);
        }
        Object mates = reader;
        List result = new ArrayList();
        boolean reverseStrand = transcript.strand() == Strand.REVERSE;
        List codingRegions =
                reverseStrand ? CollectionsKt.reversed((Iterable) codingRegions(transcript)) : codingRegions(transcript);
        for(NamedBed codingRegion : codingRegions)
        {
            Object item$iv$iv2;
            SAMRecord it;
            Iterable $receiver$iv$iv2;
            Iterable $receiver$iv2;
            Iterable bl = (Iterable) mates;
            void $i$f$mapTo = $receiver$iv2;
            Collection destination$iv$iv2 = new ArrayList();
            for(Object element$iv$iv : $receiver$iv$iv2)
            {
                it = (SAMRecord) element$iv$iv;
                boolean bl2 = false;
                if(!((long) it.getAlignmentStart() <= codingRegion.end() && (long) it.getAlignmentEnd() >= codingRegion.start()))
                {
                    continue;
                }
                destination$iv$iv2.add(element$iv$iv);
            }
            $receiver$iv2 = (List) destination$iv$iv2;
            $receiver$iv$iv2 = $receiver$iv2;
            destination$iv$iv2 = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv2, (int) 10));
            for(Object item$iv$iv2 : $receiver$iv$iv2)
            {
                it = (SAMRecord) item$iv$iv2;
                collection = destination$iv$iv2;
                boolean bl3 = false;
                object =
                        SAMCodingRecord.Companion.create$default(SAMCodingRecord.Companion, reverseStrand, (GenomeRegion) codingRegion, it, false, false, 24, null);
                collection.add(object);
            }
            $receiver$iv2 = (List) destination$iv$iv2;
            $receiver$iv$iv2 = $receiver$iv2;
            destination$iv$iv2 = new ArrayList();
            Iterable $receiver$iv$iv$iv = $receiver$iv$iv2;
            item$iv$iv2 = $receiver$iv$iv$iv.iterator();
            while(item$iv$iv2.hasNext())
            {
                NucleotideFragment nucleotideFragment;
                Object element$iv$iv$iv;
                Object element$iv$iv = element$iv$iv$iv = item$iv$iv2.next();
                SAMCodingRecord it2 = (SAMCodingRecord) element$iv$iv;
                boolean bl4 = false;
                if(mFragmentFactory.createAlignmentFragments(it2, codingRegion) == null)
                {
                    continue;
                }
                NucleotideFragment it$iv$iv = nucleotideFragment;
                destination$iv$iv2.add(it$iv$iv);
            }
            $receiver$iv2 = (List) destination$iv$iv2;
            for(Object element$iv : $receiver$iv2)
            {
                NucleotideFragment it3 = (NucleotideFragment) element$iv;
                boolean bl5 = false;
                result.add(it3);
            }
        }
        return result;

         */
    }

    private final List<SAMCodingRecord> query(
            boolean reverseStrand, GenomePosition variantRegion, NamedBed nearestCodingRegion, String bamFileName)
    {
        return Lists.newArrayList()
;
        /*
        SAMSlicer slicer = new SAMSlicer(1);
        List result = new ArrayList();
        Closeable closeable = (Closeable) samReaderFactory().open(new File(bamFileName));
        Throwable throwable = null;
        try
        {
            SamReader samReader = (SamReader) closeable;
            boolean bl = false;
            Consumer consumer = new Consumer<SAMRecord>(this, result, reverseStrand, nearestCodingRegion, slicer, variantRegion)
            {
                {
                    this$0 = sAMRecordReader;
                    $result$inlined = list;
                    $reverseStrand$inlined = bl;
                    $nearestCodingRegion$inlined = namedBed;
                    $slicer$inlined = sAMSlicer;
                    $variantRegion$inlined = genomePosition;
                }

                public final void accept(final SAMRecord samRecord)
                {
                    Intrinsics.checkParameterIsNotNull((Object) samRecord, (String) "samRecord");
                    if(SAMRecordReader.access$bothEndsInRangeOfCodingTranscripts(this$0, samRecord))
                    {
                        $result$inlined.add(SAMCodingRecord.Companion.create$default(SAMCodingRecord.Companion, $reverseStrand$inlined, (GenomeRegion) $nearestCodingRegion$inlined, samRecord, false, false, 24, null));
                    }
                    else
                    {
                        SAMRecordReader sAMRecordReader = this$0;
                        int n = SAMRecordReader.access$getAlignmentFiltered$p(sAMRecordReader);
                        SAMRecordReader.access$setAlignmentFiltered$p(sAMRecordReader, n + 1);
                    }
                }
            };
            String string = variantRegion.chromosome();
            Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "variantRegion.chromosome()");
            int n = (int) variantRegion.position();
            int n2 = (int) variantRegion.position();
            SamReader samReader2 = samReader;
            Intrinsics.checkExpressionValueIsNotNull((Object) samReader2, (String) "samReader");
            slicer.slice(string, n, n2, samReader2, consumer);
            Unit unit = Unit.INSTANCE;
        } catch (Throwable throwable2)
        {
            throwable = throwable2;
            throw throwable2;
        } finally
        {
            CloseableKt.closeFinally((Closeable) closeable, (Throwable) throwable);
        }
        return result;

         */
    }

    private final List<SAMCodingRecord> query(boolean reverseStrand, NamedBed codingRegion, String bamFileName)
    {
        return Lists.newArrayList();

        /*
        SAMSlicer slicer = new SAMSlicer(1);
        List result = new ArrayList();
        Closeable closeable = (Closeable) samReaderFactory().open(new File(bamFileName));
        Throwable throwable = null;
        try
        {
            SamReader samReader = (SamReader) closeable;
            boolean bl = false;
            Consumer consumer = new Consumer<SAMRecord>(this, result, reverseStrand, codingRegion, slicer)
            {
                {
                    this$0 = sAMRecordReader;
                    $result$inlined = list;
                    $reverseStrand$inlined = bl;
                    $codingRegion$inlined = namedBed;
                    $slicer$inlined = sAMSlicer;
                }

                public final void accept(final SAMRecord samRecord)
                {
                    Intrinsics.checkParameterIsNotNull((Object) samRecord, (String) "samRecord");
                    if(SAMRecordReader.access$bothEndsInRangeOfCodingTranscripts(this$0, samRecord))
                    {
                        $result$inlined.add(SAMCodingRecord.Companion.create$default(SAMCodingRecord.Companion, $reverseStrand$inlined, (GenomeRegion) $codingRegion$inlined, samRecord, false, false, 24, null));
                    }
                    else
                    {
                        SAMRecordReader sAMRecordReader = this$0;
                        int n = SAMRecordReader.access$getAlignmentFiltered$p(sAMRecordReader);
                        SAMRecordReader.access$setAlignmentFiltered$p(sAMRecordReader, n + 1);
                    }
                }
            };
            GenomeRegion genomeRegion = (GenomeRegion) codingRegion;
            SamReader samReader2 = samReader;
            Intrinsics.checkExpressionValueIsNotNull((Object) samReader2, (String) "samReader");
            slicer.slice(genomeRegion, samReader2, consumer);
            Unit unit = Unit.INSTANCE;
        } catch (Throwable throwable2)
        {
            throwable = throwable2;
            throw throwable2;
        } finally
        {
            CloseableKt.closeFinally((Closeable) closeable, (Throwable) throwable);
        }
        return result;

         */
    }

    private final List<NucleotideFragment> realign(NamedBed codingRegion, boolean reverseStrand, String bamFileName)
    {
        return Lists.newArrayList();

        /*
        List result = new ArrayList();
        for(SAMCodingRecord codingRecord : query(reverseStrand, codingRegion, bamFileName))
        {
            NucleotideFragment fragment;
            if(codingRecord.getIndels().contains(STOP_LOSS_ON_C))
            {
                mStopLossOnC.compute(STOP_LOSS_ON_C, realign .1.INSTANCE);
            }
            if((fragment = mFragmentFactory.createFragment(codingRecord, codingRegion)) != null)
            {
                result.add(fragment);
                continue;
            }
            for(Indel indel : codingRecord.getIndels())
            {
                if(INDEL_PON.contains(indel))
                {
                    mUnmatchedPONIndels.compute(indel, realign .2.INSTANCE);
                    continue;
                }
                mUnmatchedIndels.compute(indel, realign .3.INSTANCE);
            }
        }
        return result;

         */
    }

    private final boolean bothEndsInRangeOfCodingTranscripts(SAMRecord $receiver)
    {
        return false;

        /*
        boolean bl;
        boolean thisInRange;
        block7:
        {
            boolean bl2;
            block6:
            {
                Iterable $receiver$iv = mCodingRegions;
                if($receiver$iv instanceof Collection && ((Collection) $receiver$iv).isEmpty())
                {
                    bl2 = false;
                }
                else
                {
                    for(Object element$iv : $receiver$iv)
                    {
                        GenomeRegion it = (GenomeRegion) element$iv;
                        boolean bl3 = false;
                        if(!(Intrinsics.areEqual((Object) it.chromosome(), (Object) $receiver.getContig())
                                && (long) $receiver.getAlignmentStart() >= it.start() && (long) $receiver.getAlignmentStart() <= it.end()))
                        {
                            continue;
                        }
                        bl2 = true;
                        break block6;
                    }
                    bl2 = false;
                }
            }
            thisInRange = bl2;
            Iterable $receiver$iv = mCodingRegions;
            if($receiver$iv instanceof Collection && ((Collection) $receiver$iv).isEmpty())
            {
                bl = false;
            }
            else
            {
                for(Object element$iv : $receiver$iv)
                {
                    GenomeRegion it = (GenomeRegion) element$iv;
                    boolean bl4 = false;
                    if(!(Intrinsics.areEqual((Object) it.chromosome(), (Object) $receiver.getContig())
                            && (long) $receiver.getMateAlignmentStart() >= it.start()
                            && (long) $receiver.getMateAlignmentStart() <= it.end()))
                    {
                        continue;
                    }
                    bl = true;
                    break block7;
                }
                bl = false;
            }
        }
        boolean mateInRange = bl;
        return thisInRange && mateInRange;

         */
    }
}
