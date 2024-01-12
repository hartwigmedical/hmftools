package com.hartwig.hmftools.sage.candidate;

import static com.hartwig.hmftools.common.codon.Nucleotides.swapDnaBase;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.sage.SageConstants.MIN_CORE_DISTANCE;
import static com.hartwig.hmftools.sage.SageConstants.MIN_SOFT_CLIP_MIN_BASE_QUAL;
import static com.hartwig.hmftools.sage.evidence.RawContextCigarHandler.exceedsSoftClipLowBaseQual;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.common.ReadContext;
import com.hartwig.hmftools.sage.read.ReadContextFactory;
import com.hartwig.hmftools.sage.testutil.MutatedBases;
import com.hartwig.hmftools.sage.testutil.MutatedBasesBuilder;

import org.junit.Test;
import org.mockito.Mockito;

import htsjdk.samtools.SAMRecord;

public class RefContextConsumerTest
{
    private static final int FLANK_SIZE = 10;
    private static final ReadContextFactory READ_CONTEXT_FACTORY = new ReadContextFactory(FLANK_SIZE);
    private static final int FRAGMENT_LENGTH = 600;
    private static final int READ_LENGTH = 151;
    private static final String CHROMOSOME = CHR_1;
    private static final int RANDOM_BASES_PADDING_LENGTH = 2 * FRAGMENT_LENGTH;
    private static final String REF_BASES_STR;

    static
    {
        Random random = new Random(0);
        int middleRefBasesLength = 100;
        String initRefBases = generateRandomBases(random, RANDOM_BASES_PADDING_LENGTH);
        String middleRefBases = generateRandomBases(random, middleRefBasesLength);
        String finalRefBases = generateRandomBases(random, RANDOM_BASES_PADDING_LENGTH);
        REF_BASES_STR = initRefBases + middleRefBases + finalRefBases;
    }

    public static final IndexedBases REF_BASES = new IndexedBases(1, 0, REF_BASES_STR.getBytes());

    @Test
    public void testMnvLength()
    {
        assertEquals(1, RefContextConsumer.mnvLength(0, 0, "A".getBytes(), "CC".getBytes()));
        assertEquals(1, RefContextConsumer.mnvLength(0, 0, "AA".getBytes(), "C".getBytes()));

        assertEquals(2, RefContextConsumer.mnvLength(0, 0, "AAA".getBytes(), "CC".getBytes()));
        assertEquals(2, RefContextConsumer.mnvLength(0, 0, "AA".getBytes(), "CCC".getBytes()));
        assertEquals(2, RefContextConsumer.mnvLength(0, 0, "AAC".getBytes(), "CCC".getBytes()));

        assertEquals(3, RefContextConsumer.mnvLength(0, 0, "AAA".getBytes(), "CCC".getBytes()));
        assertEquals(3, RefContextConsumer.mnvLength(0, 0, "ACA".getBytes(), "CCCC".getBytes()));
        assertEquals(3, RefContextConsumer.mnvLength(0, 0, "ACAC".getBytes(), "CCCC".getBytes()));

        assertEquals(3, RefContextConsumer.mnvLength(0, 0, "AAAA".getBytes(), "CCCC".getBytes()));
        assertEquals(3, RefContextConsumer.mnvLength(1, 0, "AAAA".getBytes(), "CCCC".getBytes()));
        assertEquals(2, RefContextConsumer.mnvLength(2, 0, "AAAA".getBytes(), "CCCC".getBytes()));
        assertEquals(1, RefContextConsumer.mnvLength(3, 0, "AAAA".getBytes(), "CCCC".getBytes()));
    }

    @Test
    public void testLowQualSoftClips()
    {
        final byte[] baseQualities = new byte[30];

        for(int i = 0; i < baseQualities.length; ++i)
        {
            baseQualities[i] = (byte)MIN_SOFT_CLIP_MIN_BASE_QUAL;
        }

        assertFalse(exceedsSoftClipLowBaseQual(baseQualities, 0, 10));
        assertFalse(exceedsSoftClipLowBaseQual(baseQualities, 20, 10));

        for(int i = 5; i < 25; ++i)
        {
            baseQualities[i] = (byte)(MIN_SOFT_CLIP_MIN_BASE_QUAL - 1);
        }

        assertTrue(exceedsSoftClipLowBaseQual(baseQualities, 0, 10));
        assertTrue(exceedsSoftClipLowBaseQual(baseQualities, 20, 10));

        assertFalse(exceedsSoftClipLowBaseQual(baseQualities, 0, 7));
        assertFalse(exceedsSoftClipLowBaseQual(baseQualities, 23, 7));
    }

    @Test
    public void testCheckCoreExtensionAltReadPairsCheckExtended()
    {
        int startRefPos = RANDOM_BASES_PADDING_LENGTH + 1;
        int noOverlapGap = 2 * MIN_CORE_DISTANCE + 1;
        String insert = "G";
        int delLength = 1;

        CheckCoreExtensionTestCase[] testCases = new CheckCoreExtensionTestCase[] {
                new SnvInsertTestCase("NoExtensionSnvBeforeInsert", startRefPos + noOverlapGap, startRefPos, insert, false),
                new SnvInsertTestCase("NoExtensionInsertBeforeSnv", startRefPos, startRefPos + noOverlapGap, insert, false),
                new SnvDelTestCase("NoExtensionSnvBeforeDel", startRefPos + noOverlapGap + 1, startRefPos, delLength, false),
                new SnvDelTestCase("NoExtensionDelBeforeSnv", startRefPos, startRefPos + noOverlapGap + delLength, delLength, false),
                new SnvInsertTestCase("WithExtensionSnvBeforeInsert", startRefPos + 1, startRefPos, insert, true),
                new SnvInsertTestCase("WithExtensionInsertBeforeSnv", startRefPos, startRefPos + 1, insert, true),
                new SnvDelTestCase("WithExtensionSnvBeforeDel", startRefPos + 1, startRefPos, delLength, true),
                new SnvDelTestCase("WithExtensionDelBeforeSNV", startRefPos, startRefPos + 1 + delLength, delLength, true)
        };

        for(CheckCoreExtensionTestCase testCase : testCases)
            testCase.check();
    }

    private static abstract class CheckCoreExtensionTestCase
    {
        public final String Name;

        public final int IndelRefPos;
        public final int SnvRefPos;

        public final boolean ExtensionExpected;

        public CheckCoreExtensionTestCase(final String name, int indelRefPos, int snvRefPos, boolean extensionExpected)
        {
            Name = name;
            IndelRefPos = indelRefPos;
            SnvRefPos = snvRefPos;
            ExtensionExpected = extensionExpected;
        }

        protected abstract void mutate(final MutatedBasesBuilder mutatedBasesBuilder);
        protected final int indelReadIndex(int mutStartIndex)
        {
            return IndelRefPos - 1 - mutStartIndex;
        }

        private int snvReadIndex(int mutStartIndex)
        {
            int readIndex = SnvRefPos - 1 - mutStartIndex;
            if(SnvRefPos <= IndelRefPos)
            {
                return readIndex;
            }

            readIndex += indelSize();
            return readIndex;
        }

        protected abstract String indelRef();
        protected abstract String indelAlt();
        private String snvRef()
        {
            return String.valueOf(REF_BASES_STR.charAt(SnvRefPos - 1));
        }
        protected final String snvAlt()
        {
            return swapDnaBase(snvRef());
        }
        protected abstract ReadContext indelReadContext(final SAMRecord read, int mutStartIndex);
        private int indelSize()
        {
            return indelAlt().length() - indelRef().length();
        }

        private ReadContext snvReadContext(final SAMRecord read, int mutStartIndex)
        {
            return READ_CONTEXT_FACTORY.createSNVContext(SnvRefPos, snvReadIndex(mutStartIndex), read, REF_BASES);
        }

        private List<AltRead> altReadSpies()
        {
            MutatedBasesBuilder mutatedBasesBuilder = new MutatedBasesBuilder(REF_BASES_STR);
            mutate(mutatedBasesBuilder);
            MutatedBases mutatedBases = mutatedBasesBuilder.build();

            BaseRegion regionOfInterest = mutatedBases.mutIndexMutationBounds();
            int mutStartIndex = regionOfInterest.leftMidPos() - READ_LENGTH / 2;
            SAMRecord read = mutatedBases.getRead(
                    "READ_001",
                    CHROMOSOME,
                    mutStartIndex,
                    READ_LENGTH,
                    FRAGMENT_LENGTH,
                    true,
                    true);

            AltRead indelAltRead =
                    Mockito.spy(new AltRead(null, indelRef(), indelAlt(), 0, 0, false, indelReadContext(read, mutStartIndex)));
            AltRead snvAltRead = Mockito.spy(new AltRead(null, snvRef(), snvAlt(), 0, 0, false, snvReadContext(read, mutStartIndex)));

            List<AltRead> altReads = Lists.newArrayList(indelAltRead, snvAltRead);
            Collections.sort(altReads, Comparator.comparingInt(x -> x.readContext().Position));
            return altReads;
        }

        public final void check()
        {
            List<AltRead> altReads = altReadSpies();

            assertEquals(2, altReads.size());

            AltRead indelAltRead = altReads.get(0);
            AltRead snvAltRead = altReads.get(1);
            if(!indelAltRead.isIndel())
            {
                AltRead tmp = indelAltRead;
                indelAltRead = snvAltRead;
                snvAltRead = tmp;
            }

            assertTrue(indelAltRead.isIndel());
            assertFalse(snvAltRead.isIndel());

            assertFalse(indelAltRead.readContext().hasIncompleteCore());
            assertFalse(snvAltRead.readContext().hasIncompleteCore());

            BaseRegion indelCoreRegionBefore = indelAltRead.coreBaseRegion();
            BaseRegion snvCoreRegionBefore = snvAltRead.coreBaseRegion();

            assertEquals(ExtensionExpected, indelCoreRegionBefore.overlaps(snvCoreRegionBefore));

            RefContextConsumer.checkCoreExtension(altReads);

            int extendTimes = ExtensionExpected ? 1 : 0;
            Mockito.verify(indelAltRead, Mockito.times(extendTimes)).extend(Mockito.any());
            Mockito.verify(snvAltRead, Mockito.times(extendTimes)).extend(Mockito.any());
        }
    }

    private static final class SnvInsertTestCase extends CheckCoreExtensionTestCase
    {
        public final String Insert;

        public SnvInsertTestCase(final String name, int indelRefPos, int snvRefPos, final String insert, boolean extensionExpected)
        {
            super(name, indelRefPos, snvRefPos, extensionExpected);
            Insert = insert;
        }

        @Override
        protected void mutate(final MutatedBasesBuilder mutatedBasesBuilder)
        {
            mutatedBasesBuilder.insertBases(IndelRefPos, Insert);
            mutatedBasesBuilder.mutateBase(SnvRefPos, snvAlt().charAt(0));
        }

        @Override
        protected String indelRef()
        {
            return String.valueOf(REF_BASES_STR.charAt(IndelRefPos - 1));
        }
        @Override
        protected String indelAlt()
        {
            return indelRef() + Insert;
        }

        @Override
        protected ReadContext indelReadContext(final SAMRecord read, int mutStartIndex)
        {
            return READ_CONTEXT_FACTORY.createInsertContext(indelAlt(), IndelRefPos, indelReadIndex(mutStartIndex), read.getReadBases(), REF_BASES);
        }

        @Override
        public String toString()
        {
            return "SnvInsertTestCase{" +
                    "Name='" + Name + '\'' +
                    ", IndelRefPos=" + IndelRefPos +
                    ", SnvRefPos=" + SnvRefPos +
                    ", Insert='" + Insert + '\'' +
                    ", ExtensionExpected=" + ExtensionExpected +
                    '}';
        }
    }

    private static final class SnvDelTestCase extends CheckCoreExtensionTestCase
    {
        public final int DelLength;

        public SnvDelTestCase(final String name, int indelRefPos, int snvRefPos, final int delLength, boolean extensionExpected)
        {
            super(name, indelRefPos, snvRefPos, extensionExpected);
            DelLength = delLength;
        }

        @Override
        protected void mutate(final MutatedBasesBuilder mutatedBasesBuilder)
        {
            mutatedBasesBuilder.delBases(IndelRefPos + 1, DelLength);
            mutatedBasesBuilder.mutateBase(SnvRefPos, snvAlt().charAt(0));
        }

        @Override
        protected String indelRef()
        {
            return REF_BASES_STR.substring(IndelRefPos - 1, IndelRefPos + DelLength);
        }
        @Override
        protected String indelAlt()
        {
            return String.valueOf(REF_BASES_STR.charAt(IndelRefPos - 1));
        }

        @Override
        protected ReadContext indelReadContext(final SAMRecord read, int mutStartIndex)
        {
            return READ_CONTEXT_FACTORY.createDelContext(indelRef(), IndelRefPos, indelReadIndex(mutStartIndex), read.getReadBases(), REF_BASES);
        }

        @Override
        public String toString()
        {
            return "SnvDelTestCase{" +
                    "Name='" + Name + '\'' +
                    ", IndelRefPos=" + IndelRefPos +
                    ", SnvRefPos=" + SnvRefPos +
                    ", DelLength=" + DelLength +
                    ", ExtensionExpected=" + ExtensionExpected +
                    '}';
        }
    }
}
