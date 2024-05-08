package com.hartwig.hmftools.esvee.depth;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.SvVcfTags.TOTAL_FRAGS;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.SamRecordTestUtils.setReadFlag;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.test.ReadIdGenerator;

import org.junit.Assert;
import org.junit.Test;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

public class DepthAnnotationTest
{
    private final DepthTask mDepthTask;
    private int mNextVariantId;
    private final ReadIdGenerator mReadIdGen;

    private static final String TEST_SAMPLE_ID = "SAMPLE";

    public DepthAnnotationTest()
    {
        mNextVariantId = 0;
        mReadIdGen = new ReadIdGenerator();

        Map<String,Integer> sampleVcfGenotypeIds = Maps.newHashMap();
        sampleVcfGenotypeIds.put(TEST_SAMPLE_ID, 0);
        DepthConfig config = new DepthConfig(0.1, 1000);
        config.Samples.add(TEST_SAMPLE_ID);
        mDepthTask = new DepthTask(CHR_1, config, sampleVcfGenotypeIds);
    }

    @Test
    public void testReferenceDepth()
    {
        mDepthTask.reset();

        VariantContext var1 = createVariantContext(nextVariantId(), 1000);
        VariantContext var2 = createVariantContext(nextVariantId(), 1100);
        VariantContext var3 = createVariantContext(nextVariantId(), 1200);

        List<VariantContext> variants = Lists.newArrayList(var1, var2, var3);

        mDepthTask.addVariants(variants);

        mDepthTask.addSliceVariants(mDepthTask.variantInfos());

        // irrelevant
        SAMRecord read1 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 900, CHR_2, 100, "100M", true, false,
                false, "");

        SAMRecord read2 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 1250, CHR_2, 100, "100M", true, false,
                false, "");

        mDepthTask.processSamRecord(read1);
        mDepthTask.processSamRecord(read2);
        assertTrue(mDepthTask.readGroups().isEmpty());

        checkRefSupport(1000, 0, 0);

        // reads processed immediately for overlapping overage
        SAMRecord read3 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 950, CHR_2, 100, "100M", true, false,
                false, "");

        mDepthTask.processSamRecord(read3);

        assertTrue(mDepthTask.readGroups().isEmpty());

        checkRefSupport(1000, 1, 0);
        checkRefSupport(1100, 0, 0);

        // a read straddling all three variants
        SAMRecord read4 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 800, CHR_1, 1300, "100M", true, false,
                false, "");

        SAMRecord mate4 = createSamRecord(
                read4.getReadName(), CHR_1, 1300, CHR_1, 800, "100M", false, true,
                false, "");
        mate4.setReadNegativeStrandFlag(true);

        mDepthTask.processSamRecord(read4);

        Assert.assertEquals(1, mDepthTask.readGroups().size());

        mDepthTask.processSamRecord(mate4);
        assertTrue(mDepthTask.readGroups().isEmpty());

        checkRefSupport(1000, 1, 1);
        checkRefSupport(1100, 0, 1);
        checkRefSupport(1200, 0, 1);

        // confirm that junction-matching reads are ignored
        SAMRecord read5 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 1100, CHR_2, 100, "5S145M", true, false,
                false, "");

        mDepthTask.processSamRecord(read5);

        checkRefSupport(1100, 1, 1);
        checkRefSupport(1200, 1, 1);

        SAMRecord read6 = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 960, CHR_2, 100, "141M10S", true, false,
            false, "");

        mDepthTask.processSamRecord(read6);

        checkRefSupport(1000, 2, 1);
        checkRefSupport(1100, 1, 1);
        checkRefSupport(1200, 1, 1);

        // fragments with local supplementaries wait until all are received

        // test max fragment caps
        for(int i = 0; i < 15; ++i)
        {
            SAMRecord read = createSamRecord(
                    mReadIdGen.nextId(), CHR_1, 1150, CHR_2, 100, "100M", true, false,
                    false, "");

            mDepthTask.processSamRecord(read);
        }

        checkRefSupport(1200, 9, 1);
        Assert.assertEquals(2, mDepthTask.sliceRegionState().UncappedVariants.size());
    }

    @Test
    public void testPositionProgression()
    {
        mDepthTask.reset();

        int startPosition = 1000;
        int posGap = 200;
        List<VariantContext> variants = Lists.newArrayList();

        for(int i = 0; i < 20; ++i)
        {
            VariantContext var = createVariantContext(nextVariantId(), startPosition + i * posGap);
            variants.add(var);
        }

        mDepthTask.addVariants(variants);
        mDepthTask.addSliceVariants(mDepthTask.variantInfos());

        // check slice state as reads are received
        SAMRecord read = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 950, CHR_2, 100, "100M", true, false,
                false, "");

        mDepthTask.processSamRecord(read);

        Assert.assertEquals(variants.size(), mDepthTask.sliceRegionState().UncappedVariants.size());
        Assert.assertEquals(0, mDepthTask.sliceRegionState().MinPositionIndex);

        read = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 1150, CHR_2, 100, "100M", true, false,
                false, "");

        mDepthTask.processSamRecord(read);

        Assert.assertEquals(1, mDepthTask.sliceRegionState().MinPositionIndex);

        // process ahead and then cap out some variants
        read = createSamRecord(
                mReadIdGen.nextId(), CHR_1, 2150, CHR_2, 100, "100M", true, false,
                false, "");

        mDepthTask.processSamRecord(read);

        Assert.assertEquals(6, mDepthTask.sliceRegionState().MinPositionIndex);

        for(int i = 0; i < 15; ++i)
        {
            SAMRecord read2 = createSamRecord(
                    mReadIdGen.nextId(), CHR_1, 2150, CHR_2, 100, "100M", true, false,
                    false, "");

            mDepthTask.processSamRecord(read2);
        }

        Assert.assertEquals(5, mDepthTask.sliceRegionState().MinPositionIndex);
        Assert.assertEquals(variants.size() - 1, mDepthTask.sliceRegionState().UncappedVariants.size());
    }

    private void checkRefSupport(int varPosition, int refSupport, int refPairSupport)
    {
        VariantInfo variant = mDepthTask.variantInfos().stream().filter(x -> x.Position == varPosition).findFirst().orElse(null);
        assertNotNull(variant);

        Assert.assertEquals(refSupport, variant.SampleSupportCounts[0].RefSupport);
        Assert.assertEquals(refPairSupport, variant.SampleSupportCounts[0].RefPairSupport);
    }

    private String nextVariantId() { return format("%03d", mNextVariantId++); }

    private static SAMRecord createSamRecord(
            final String readId, final String chromosome, int readStart, final String mateChr, int mateStart,
            final String cigar, boolean firstInPair, boolean isReversed, boolean isSupplementary, final String suppData)
    {
        SAMRecordSetBuilder recordBuilder = new SAMRecordSetBuilder();
        recordBuilder.setUnmappedHasBasesAndQualities(false);

        SAMRecord record = recordBuilder.addFrag(
                readId, 1, readStart, isReversed, false, cigar, "", 10, false);

        /*
        record.setReadBases(readBases.getBytes());

        final byte[] qualities = new byte[readBases.length()];

        for(int i = 0; i < readBases.length(); ++i)
            qualities[i] = (byte)baseQual;

        record.setBaseQualities(qualities);
        */

        record.setReferenceName(chromosome);
        record.setReferenceIndex(chromosomeToIndex(chromosome));

        record.setMateReferenceName(mateChr);
        record.setMateAlignmentStart(mateStart);
        record.setMateReferenceIndex(chromosomeToIndex(mateChr));

        int flags = 0;

        flags = setReadFlag(flags, SAMFlag.READ_PAIRED);
        flags = setReadFlag(flags, SAMFlag.PROPER_PAIR);

        if(isReversed)
            flags = setReadFlag(flags, SAMFlag.READ_REVERSE_STRAND);

        if(firstInPair)
            flags = setReadFlag(flags, SAMFlag.FIRST_OF_PAIR);
        else
            flags = setReadFlag(flags, SAMFlag.SECOND_OF_PAIR);

        // if(secondary)
        //    flags = setReadFlag(flags, SAMFlag.SECONDARY_ALIGNMENT);

        if(isSupplementary)
            flags = setReadFlag(flags, SAMFlag.SUPPLEMENTARY_ALIGNMENT);

        record.setFlags(flags);
        return record;
    }

    private static int chromosomeToIndex(final String chromosome)
    {
        return Integer.parseInt(chromosome) - 1;
    }

    public static VariantContext createVariantContext(final String vcfId, int position)
    {
        // required - VF, VARIANT_FRAGMENT_BREAKEND_COVERAGE
        VariantContextBuilder builder = new VariantContextBuilder();

        List<Allele> alleles = Lists.newArrayList();

        alleles.add(Allele.create("A", true));
        alleles.add(Allele.create("T[2:500[", false));

        Map<String,Object> commonAttributes = Maps.newHashMap();
        commonAttributes.put(TOTAL_FRAGS, 1);

        Map<String,Object> genotypeAttributes = Maps.newHashMap();

        // defaults to indicate a somatic variant
        genotypeAttributes.put(TOTAL_FRAGS, 1);

        Genotype genotype = new GenotypeBuilder()
                .attributes(genotypeAttributes)
                .name(TEST_SAMPLE_ID)
                .DP(-1)
                .noAD()
                .noPL()
                .GQ(-1)
                .make();

        GenotypesContext genotypesContext = GenotypesContext.create(genotype);

        double logError = -(37 / 10.0);
        return builder
                .source("SOURCE")
                .id(vcfId)
                .chr(CHR_1)
                .start(position)
                .stop(position)
                .alleles(alleles)
                .genotypes(genotypesContext)
                .attributes(commonAttributes)
                .log10PError(logError)
                .unfiltered()
                .make(true);
    }





}
