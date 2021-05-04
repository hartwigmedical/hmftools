package com.hartwig.hmftools.lilac.sam;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.samtools.CigarHandler;
import com.hartwig.hmftools.common.samtools.CigarTraversal;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

public class SAMCodingRecord
{
    private final String id;
    private final int softClippedStart;
    private final int softClippedEnd;
    private final List<Indel> indels;
    private final int positionStart;
    private final int positionEnd;
    private final int readStart;
    private final int readEnd;
    private final SAMRecord record;
    private final boolean reverseStrand;

    public SAMCodingRecord(final String id, int softClippedStart, int softClippedEnd, final List<Indel> indels2, int positionStart,
            int positionEnd, int readStart, int readEnd, final SAMRecord record, boolean reverseStrand)
    {
        this.id = id;
        this.softClippedStart = softClippedStart;
        this.softClippedEnd = softClippedEnd;
        this.indels = indels2;
        this.positionStart = positionStart;
        this.positionEnd = positionEnd;
        this.readStart = readStart;
        this.readEnd = readEnd;
        this.record = record;
        this.reverseStrand = reverseStrand;
    }

    public final int maxIndelSize()
    {
        return 1;
        // TODO
        
        /*
        void var3_3;
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        Iterable iterable = $receiver$iv = (Iterable) this.indels;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            void it;
            Indel indel = (Indel) item$iv$iv;
            Collection collection = destination$iv$iv;
            boolean bl = false;
            int n = it.getLength();
            Integer n2 = Math.abs(n);
            collection.add(n2);
        }
        Integer n = (Integer) CollectionsKt.max((Iterable) ((List) var3_3));
        return n != null ? n : 0;
        
         */
    }

    public final boolean containsSoftClip()
    {
        return this.softClippedStart > 0 || this.softClippedEnd > 0;
    }

    public final boolean containsIndel()
    {
        Collection collection = this.indels;
        return !collection.isEmpty();
    }

    public final char[] codingRegionRead(boolean reverseCompliment)
    {
        return null;
        
        /*
        char[] cArray;
        if(reverseCompliment)
        {
            void $receiver$iv$iv;
            char[] $receiver$iv;
            char[] cArray2 = $receiver$iv = this.forwardRead();
            Collection destination$iv$iv = new ArrayList($receiver$iv.length);
            int n = ((void) $receiver$iv$iv).length;
            for(int i = 0; i < n; ++i)
            {
                void it;
                void item$iv$iv;
                void var8_8 = item$iv$iv = $receiver$iv$iv[i];
                Collection collection = destination$iv$iv;
                boolean bl = false;
                Character c = Character.valueOf(Companion.reverseCompliment((char) it));
                collection.add(c);
            }
            cArray = CollectionsKt.toCharArray((Collection) CollectionsKt.reversed((Iterable) ((List) destination$iv$iv)));
        }
        else
        {
            cArray = this.forwardRead();
        }
        return cArray;
        
         */
    }

    public final int[] codingRegionQuality(boolean reverseCompliment)
    {
        return null;

//        return reverseCompliment
//                ? CollectionsKt.toIntArray((Collection) ArraysKt.reversed((int[]) this.forwardQuality()))
//                : this.forwardQuality();
    }

    private final char[] forwardRead()
    {
        return null;
        
        /*
        void var3_4;
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        int n = this.readStart;
        Iterable iterable = $receiver$iv = (Iterable) new IntRange(n, this.readEnd);
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        Iterator iterator = $receiver$iv$iv.iterator();
        while(iterator.hasNext())
        {
            void it;
            int item$iv$iv;
            int n2 = item$iv$iv = ((IntIterator) iterator).nextInt();
            Collection collection = destination$iv$iv;
            boolean bl = false;
            Character c = Character.valueOf((char) this.record.getReadBases()[it]);
            collection.add(c);
        }
        return CollectionsKt.toCharArray((Collection) ((List) var3_4));
        
         */
    }

    private final int[] forwardQuality()
    {
        return null;
        
        /*
        void var3_4;
        void $receiver$iv$iv;
        Iterable $receiver$iv;
        int n = this.readStart;
        Iterable iterable = $receiver$iv = (Iterable) new IntRange(n, this.readEnd);
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        Iterator iterator = $receiver$iv$iv.iterator();
        while(iterator.hasNext())
        {
            void it;
            int item$iv$iv;
            int n2 = item$iv$iv = ((IntIterator) iterator).nextInt();
            Collection collection = destination$iv$iv;
            boolean bl = false;
            Integer n3 = this.record.getBaseQualities()[it];
            collection.add(n3);
        }
        return CollectionsKt.toIntArray((Collection) ((List) var3_4));
        
         */
    }

    public final List<SAMCodingRecord> alignmentsOnly()
    {
        return Lists.newArrayList();

        /*
        Object object;
        GenomeRegion it;
        Collection collection;
        Iterable $receiver$iv$iv;
        GenomeRegion outerRegion =
                GenomeRegions.create((String) this.record.getContig(), (long) this.positionStart, (long) this.positionEnd);
        List list = this.record.getAlignmentBlocks();
        Intrinsics.checkExpressionValueIsNotNull((Object) list, (String) "record.alignmentBlocks");
        Iterable $receiver$iv = list;
        Iterable iterable = $receiver$iv;
        Collection destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            AlignmentBlock alignmentBlock = (AlignmentBlock) item$iv$iv;
            collection = destination$iv$iv;
            boolean bl = false;
            String string = this.record.getContig();
            void v2 = it;
            Intrinsics.checkExpressionValueIsNotNull((Object) v2, (String) "it");
            object = GenomeRegions.create((String) string, (long) v2.getReferenceStart(), (long) (
                    (long) it.getReferenceStart() + (long) it.getLength() - 1L));
            collection.add(object);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList();
        for(Object element$iv$iv : $receiver$iv$iv)
        {
            it = (GenomeRegion) element$iv$iv;
            boolean bl = false;
            if(!outerRegion.overlaps(it))
            {
                continue;
            }
            destination$iv$iv.add(element$iv$iv);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            it = (GenomeRegion) item$iv$iv;
            collection = destination$iv$iv;
            boolean bl = false;
            long $i$f$filterTo = it.start();
            long l = outerRegion.start();
            long l2 = Math.max($i$f$filterTo, l);
            $i$f$filterTo = it.end();
            l = outerRegion.end();
            object = GenomeRegions.create((String) this.record.getContig(), (long) l2, (long) Math.min($i$f$filterTo, l));
            collection.add(object);
        }
        $receiver$iv = (List) destination$iv$iv;
        $receiver$iv$iv = $receiver$iv;
        destination$iv$iv = new ArrayList(CollectionsKt.collectionSizeOrDefault((Iterable) $receiver$iv, (int) 10));
        for(Object item$iv$iv : $receiver$iv$iv)
        {
            it = (GenomeRegion) item$iv$iv;
            collection = destination$iv$iv;
            boolean bl = false;
            GenomeRegion genomeRegion = it;
            Intrinsics.checkExpressionValueIsNotNull((Object) genomeRegion, (String) "it");
            object = Companion.create(this.reverseStrand, genomeRegion, this.record, false, false);
            collection.add(object);
        }
        List result = (List) destination$iv$iv;
        return result;
        
         */
    }

    public final String getId()
    {
        return this.id;
    }

    public final int getSoftClippedStart()
    {
        return this.softClippedStart;
    }

    public final int getSoftClippedEnd()
    {
        return this.softClippedEnd;
    }

    public final List<Indel> getIndels()
    {
        return this.indels;
    }

    public final int getPositionStart()
    {
        return this.positionStart;
    }

    public final int getPositionEnd()
    {
        return this.positionEnd;
    }

    public final int getReadStart()
    {
        return this.readStart;
    }

    public final int getReadEnd()
    {
        return this.readEnd;
    }

    public final SAMRecord getRecord()
    {
        return this.record;
    }

    public final boolean getReverseStrand()
    {
        return this.reverseStrand;
    }

    public static SAMCodingRecord create(boolean reverseStrand, final GenomeRegion codingRegion, final SAMRecord record,
            boolean includeSoftClips, boolean includeIndels)
    {
        return null;

        /*
        int n;
        Intrinsics.checkParameterIsNotNull((Object) codingRegion, (String) "codingRegion");
        Intrinsics.checkParameterIsNotNull((Object) record, (String) "record");
        int softClipStart = this.softClipStart(record);
        int softClipEnd = this.softClipEnd(record);
        int alignmentStart = record.getAlignmentStart();
        int alignmentEnd = record.getAlignmentEnd();
        int recordStart = alignmentStart - softClipStart;
        int recordEnd = alignmentEnd + softClipEnd;
        int n2 = (int) codingRegion.start();
        int positionStart = Math.max(n2, alignmentStart);
        int n3 = (int) codingRegion.end();
        int positionEnd = Math.min(n3, alignmentEnd);
        int readIndexStart = record.getReadPositionAtReferencePosition(positionStart, true) - 1;
        int readIndexEnd = record.getReadPositionAtReferencePosition(positionEnd, true) - 1;
        positionStart = record.getReferencePositionAtReadPosition(readIndexStart + 1);
        positionEnd = record.getReferencePositionAtReadPosition(readIndexEnd + 1);
        if(positionStart == alignmentStart && softClipStart > 0 && includeSoftClips)
        {
            n = (int) codingRegion.start();
            int earliestStart = Math.max(n, recordStart);
            readIndexStart = readIndexStart - positionStart + earliestStart;
            positionStart = earliestStart;
        }
        if(positionEnd == alignmentEnd && softClipEnd > 0 && includeSoftClips)
        {
            n = (int) codingRegion.end();
            int latestEnd = Math.min(n, recordEnd);
            readIndexEnd = readIndexEnd + latestEnd - positionEnd;
            positionEnd = latestEnd;
        }
        n = alignmentStart - positionStart;
        int n4 = 0;
        int softClippedStart = Math.max(n, n4);
        n4 = 0;
        int n5 = positionEnd - alignmentEnd;
        int softClippedEnd = Math.max(n4, n5);
        List<Indel> indels2 = includeIndels ? this.indels(positionStart, positionEnd, record) : CollectionsKt.emptyList();
        String string = record.getReadName();
        Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "record.readName");
        return new SAMCodingRecord(string, softClippedStart, softClippedEnd, indels2, positionStart, positionEnd, readIndexStart, readIndexEnd, record, reverseStrand);

         */
    }

    private static List<Indel> indels(int startPosition, int endPosition, SAMRecord record)
    {
        return Lists.newArrayList();

        /*
        List indels2 = new ArrayList();
        CigarHandler handler2 = new CigarHandler(startPosition, endPosition, indels2)
        {
            public void handleInsert(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
            {
                Intrinsics.checkParameterIsNotNull((Object) record, (String) "record");
                Intrinsics.checkParameterIsNotNull((Object) element, (String) "element");
                int n = refPosition;
                if(this.$startPosition <= n && this.$endPosition >= n)
                {
                    char base = (char) record.getReadBases()[readIndex];
                    String string = record.getContig();
                    Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "record.contig");
                    String string2 = String.valueOf(base);
                    String string3 = record.getReadString();
                    Intrinsics.checkExpressionValueIsNotNull((Object) string3, (String) "record.readString");
                    String string4 = string3;
                    int n2 = readIndex + element.getLength() + 1;
                    String string5 = string4;
                    if(string5 == null)
                    {
                        throw new TypeCastException("null cannot be cast to non-null type java.lang.String");
                    }
                    String string6 = string5.substring(readIndex, n2);
                    Intrinsics.checkExpressionValueIsNotNull((Object) string6, (String) "(this as java.lang.Strin\u2026ing(startIndex, endIndex)");
                    Indel
                            insert = new Indel(string, refPosition, string2, string6);
                    this.$indels.add(insert);
                }
            }

            public void handleDelete(final SAMRecord record, final CigarElement element, int readIndex, int refPosition)
            {
                Intrinsics.checkParameterIsNotNull((Object) record, (String) "record");
                Intrinsics.checkParameterIsNotNull((Object) element, (String) "element");
                int n = refPosition;
                if(this.$startPosition <= n && this.$endPosition >= n)
                {
                    char base = (char) record.getReadBases()[readIndex];
                    String string = record.getContig();
                    Intrinsics.checkExpressionValueIsNotNull((Object) string, (String) "record.contig");
                    char c = base;
                    String string2 = StringsKt.repeat((CharSequence) "N", (int) element.getLength());
                    Indel delete = new Indel(string, refPosition, String.valueOf(c) + string2, String.valueOf(base));
                    this.$indels.add(delete);
                }
            }

            {
                this.$startPosition = $captured_local_variable$0;
                this.$endPosition = $captured_local_variable$1;
                this.$indels = $captured_local_variable$2;
            }
        };
        CigarTraversal.traverseCigar((SAMRecord) record, (CigarHandler) handler2);
        return indels2;

         */
    }

    private static int softClipStart(SAMRecord $receiver)
    {
        return 1;

        /*
        int n;
        Cigar cigar = $receiver.getCigar();
        Intrinsics.checkExpressionValueIsNotNull((Object) cigar, (String) "this.cigar");
        CigarElement cigarElement = cigar.getFirstCigarElement();
        Intrinsics.checkExpressionValueIsNotNull((Object) cigarElement, (String) "this.cigar.firstCigarElement");
        if(cigarElement.getOperator() == CigarOperator.S)
        {
            Cigar cigar2 = $receiver.getCigar();
            Intrinsics.checkExpressionValueIsNotNull((Object) cigar2, (String) "this.cigar");
            CigarElement cigarElement2 = cigar2.getFirstCigarElement();
            Intrinsics.checkExpressionValueIsNotNull((Object) cigarElement2, (String) "this.cigar.firstCigarElement");
            n = cigarElement2.getLength();
        }
        else
        {
            n = 0;
        }
        return n;

         */
    }

    private final int softClipEnd(SAMRecord $receiver)
    {
        return 1;

        /*
        int n;
        Cigar cigar = $receiver.getCigar();
        Intrinsics.checkExpressionValueIsNotNull((Object) cigar, (String) "this.cigar");
        CigarElement cigarElement = cigar.getLastCigarElement();
        Intrinsics.checkExpressionValueIsNotNull((Object) cigarElement, (String) "this.cigar.lastCigarElement");
        if(cigarElement.getOperator() == CigarOperator.S)
        {
            Cigar cigar2 = $receiver.getCigar();
            Intrinsics.checkExpressionValueIsNotNull((Object) cigar2, (String) "this.cigar");
            CigarElement cigarElement2 = cigar2.getLastCigarElement();
            Intrinsics.checkExpressionValueIsNotNull((Object) cigarElement2, (String) "this.cigar.lastCigarElement");
            n = cigarElement2.getLength();
        }
        else
        {
            n = 0;
        }
        return n;

         */
    }

    public final char reverseCompliment(char $receiver)
    {
        switch($receiver)
        {
            case 'G':
            {
                return 'C';
            }
            case 'A':
            {
                return 'T';
            }
            case 'T':
            {
                return 'A';
            }
            case 'C':
            {
                return 'G';
            }
        }
        return $receiver;
    }
}
