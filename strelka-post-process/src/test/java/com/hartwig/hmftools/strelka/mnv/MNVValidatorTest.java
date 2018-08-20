package com.hartwig.hmftools.strelka.mnv;

import static com.hartwig.hmftools.strelka.StrelkaPostProcessApplication.generateOutputHeader;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.collect.Streams;
import com.google.common.io.Resources;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

public class MNVValidatorTest {
    private static final Logger LOGGER = LogManager.getLogger(MNVValidatorTest.class);

    private static final int GAP_SIZE = 1;
    private static final File VCF_FILE = new File(Resources.getResource("mnvs.vcf").getPath());
    private static final VCFFileReader VCF_FILE_READER = new VCFFileReader(VCF_FILE, false);
    private static final VCFHeader VCF_OUTPUT_HEADER = generateOutputHeader(VCF_FILE_READER.getFileHeader(), "TUMOR");
    private static final MNVMerger MNV_MERGER = ImmutableMNVMerger.of(VCF_OUTPUT_HEADER);
    private static final List<VariantContext> VARIANTS = Streams.stream(VCF_FILE_READER).collect(Collectors.toList());

    // MIVO: variants at positions: 1  3(alts: A,T)  =>  possible mnvs: (1,3A) (1,3T) => actual mnv: (1,3A) => non-mnv variant: 3T
    @Test
    public void correctlyDetectsNonMnvVariant() {
        final PotentialMNVRegion region =
                PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(), Lists.newArrayList(VARIANTS.get(7), VARIANTS.get(9)), GAP_SIZE);
        assertEquals(2, region.potentialMnvs().size());
        final List<VariantContext> nonMnvs = MNVRegionValidator.nonMnvVariants(region, Sets.newHashSet(region.potentialMnvs().get(0)));
        assertEquals(1, nonMnvs.size());
        assertEquals("G", nonMnvs.get(0).getReference().getBaseString());
        assertEquals(1, nonMnvs.get(0).getAlternateAlleles().size());
        assertEquals("T", nonMnvs.get(0).getAlternateAllele(0).getBaseString());
    }

    // MIVO: variants at positions: 1  3(alts: A,T)  =>  possible mnvs: (1,3A) (1,3T) => actual mnv: none => non-mnv variants: 1  3(alts: A,T)
    @Test
    public void correctlyDetectsMultiAltNonMnvVariant() {
        final PotentialMNVRegion region =
                PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(), Lists.newArrayList(VARIANTS.get(7), VARIANTS.get(9)), GAP_SIZE);
        assertEquals(2, region.potentialMnvs().size());
        final List<VariantContext> nonMnvs = MNVRegionValidator.nonMnvVariants(region, Sets.newHashSet());
        assertEquals(2, nonMnvs.size());
        assertEquals("C", nonMnvs.get(0).getReference().getBaseString());
        assertEquals(1, nonMnvs.get(0).getAlternateAlleles().size());
        assertEquals("T", nonMnvs.get(0).getAlternateAllele(0).getBaseString());
        assertEquals("G", nonMnvs.get(1).getReference().getBaseString());
        assertEquals(2, nonMnvs.get(1).getAlternateAlleles().size());
        assertEquals("A", nonMnvs.get(1).getAlternateAllele(0).getBaseString());
        assertEquals("T", nonMnvs.get(1).getAlternateAllele(1).getBaseString());
    }

    // MIVO: variants: 170756001 (C->T),  170756002 (G->T)
    // MIVO:    2 records contain both;
    @Test
    public void correctOutputForMnvOf2InRegionOf2() {
        final PotentialMNVRegion region = PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(),
                Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12)),
                GAP_SIZE);
        final SAMRecord recordWith1and2 = TestUtils.buildSamRecord(170755999, "6M", "CATTAG");
        final List<SAMRecord> records = Lists.newArrayList(recordWith1and2, recordWith1and2);
        final MNVRegionValidator validator = MNVValidator.validateMNVs(records.iterator(), region);
        final List<VariantContext> outputVariants = MNVValidator.outputVariants(validator, MNV_MERGER);
        LOGGER.info(outputVariants);
        assertEquals(1, outputVariants.size());
        assertEquals("CG", outputVariants.get(0).getReference().getBaseString());
        assertEquals(1, outputVariants.get(0).getAlternateAlleles().size());
        assertEquals("TT", outputVariants.get(0).getAlternateAllele(0).getBaseString());
    }

    // MIVO: variants: 170756001 (C->T),  170756002 (G->T)
    // MIVO:    1 record contains first; 1 record contains second, 1 record contains both
    @Test
    public void correctOutputForNoMnvOfInRegionOf2() {
        final PotentialMNVRegion region = PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(),
                Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12)),
                GAP_SIZE);
        final SAMRecord recordWith1and2 = TestUtils.buildSamRecord(170755999, "6M", "CATTAG");
        final SAMRecord recordWith1 = TestUtils.buildSamRecord(170755999, "6M", "CATGAG");
        final SAMRecord recordWith2 = TestUtils.buildSamRecord(170755999, "6M", "CACTAG");
        final List<SAMRecord> records = Lists.newArrayList(recordWith1, recordWith2, recordWith1and2);
        final MNVRegionValidator validator = MNVValidator.validateMNVs(records.iterator(), region);
        final List<VariantContext> outputVariants = MNVValidator.outputVariants(validator, MNV_MERGER);
        LOGGER.info(outputVariants);
        assertEquals(2, outputVariants.size());
        assertEquals(VARIANTS.get(11), outputVariants.get(0));
        assertEquals(VARIANTS.get(12), outputVariants.get(1));
    }

    // MIVO: variants: 170756001 (C->T),  170756002 (G->T),  170756003 (A->T)
    // MIVO:    2 records contain all 3
    @Test
    public void correctOutputForMnvOf3InRegionOf3() {
        final PotentialMNVRegion region = PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(),
                Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12), VARIANTS.get(13)),
                GAP_SIZE);
        final SAMRecord recordWith1and2and3 = TestUtils.buildSamRecord(170756000, "6M", "ATTTGC");
        final List<SAMRecord> records = Lists.newArrayList(recordWith1and2and3, recordWith1and2and3);
        final MNVRegionValidator validator = MNVValidator.validateMNVs(records.iterator(), region);
        final List<VariantContext> outputVariants = MNVValidator.outputVariants(validator, MNV_MERGER);
        assertEquals(1, outputVariants.size());
        assertEquals("CGA", outputVariants.get(0).getReference().getBaseString());
        assertEquals(1, outputVariants.get(0).getAlternateAlleles().size());
        assertEquals("TTT", outputVariants.get(0).getAlternateAllele(0).getBaseString());
    }

    // MIVO: variants: 170756001 (C->T),  170756002 (G->T),  170756003 (A->T)
    // MIVO:    1st record contains all 3; 2nd record contains 1 and 2
    @Test
    public void correctOutputForMnvOf2InRegionOf3() {
        final PotentialMNVRegion region = PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(),
                Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12), VARIANTS.get(13)),
                GAP_SIZE);
        final SAMRecord recordWith1and2 = TestUtils.buildSamRecord(170755999, "6M", "CATTAG");
        final SAMRecord recordWith1and2and3 = TestUtils.buildSamRecord(170756000, "6M", "ATTTGC");
        final List<SAMRecord> records = Lists.newArrayList(recordWith1and2, recordWith1and2and3);
        final MNVRegionValidator validator = MNVValidator.validateMNVs(records.iterator(), region);
        final List<VariantContext> outputVariants = MNVValidator.outputVariants(validator, MNV_MERGER);
        assertEquals(2, outputVariants.size());
        assertEquals("CG", outputVariants.get(0).getReference().getBaseString());
        assertEquals(1, outputVariants.get(0).getAlternateAlleles().size());
        assertEquals("TT", outputVariants.get(0).getAlternateAllele(0).getBaseString());
        assertEquals(VARIANTS.get(13), outputVariants.get(1));
    }

    // MIVO: variants: 15 (C->T),  15 (CA->C),  17 (G->A)
    // MIVO:    5 records contain 1 and 3; 1 record contains 2 and 3
    @Test
    public void correctOutputForMnvOf2InRegionOf3With2AtSamePos() {
        final PotentialMNVRegion region = PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(),
                Lists.newArrayList(VARIANTS.get(15), VARIANTS.get(16), VARIANTS.get(17)),
                GAP_SIZE);
        final SAMRecord recordWith1and3 = TestUtils.buildSamRecord(14, "6M", "CTCAGC");
        final SAMRecord recordWith2and3 = TestUtils.buildSamRecord(14, "2M1D3M", "CCAGC");
        final List<SAMRecord> records =
                Lists.newArrayList(recordWith1and3, recordWith1and3, recordWith1and3, recordWith1and3, recordWith1and3, recordWith2and3);
        final MNVRegionValidator validator = MNVValidator.validateMNVs(records.iterator(), region);
        final List<VariantContext> outputVariants = MNVValidator.outputVariants(validator, MNV_MERGER);
        assertEquals(2, outputVariants.size());
        assertEquals(VARIANTS.get(16), outputVariants.get(0));
        assertEquals("CCG", outputVariants.get(1).getReference().getBaseString());
        assertEquals(1, outputVariants.get(1).getAlternateAlleles().size());
        assertEquals("TCA", outputVariants.get(1).getAlternateAllele(0).getBaseString());
    }

    // MIVO: variants: 170755901 (C->T),  170755903 (G->A,T)
    // MIVO:    5 records contain 1 and (G->A); 1 record contains 1 and (G->T)
    @Test
    public void correctOutputForMnvOf2With2Alts() {
        final PotentialMNVRegion region =
                PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(), Lists.newArrayList(VARIANTS.get(7), VARIANTS.get(9)), GAP_SIZE);
        final SAMRecord recordWith1and2 = TestUtils.buildSamRecord(170755900, "6M", "CTCAAG");
        final SAMRecord recordWith1and3 = TestUtils.buildSamRecord(170755900, "6M", "CTCTAG");
        final List<SAMRecord> records =
                Lists.newArrayList(recordWith1and2, recordWith1and2, recordWith1and2, recordWith1and2, recordWith1and2, recordWith1and3);
        final MNVRegionValidator validator = MNVValidator.validateMNVs(records.iterator(), region);
        final List<VariantContext> outputVariants = MNVValidator.outputVariants(validator, MNV_MERGER);
        assertEquals(2, outputVariants.size());
        assertEquals("CCG", outputVariants.get(0).getReference().getBaseString());
        assertEquals(1, outputVariants.get(0).getAlternateAlleles().size());
        assertEquals("TCA", outputVariants.get(0).getAlternateAllele(0).getBaseString());
        assertEquals("G", outputVariants.get(1).getReference().getBaseString());
        assertEquals(1, outputVariants.get(1).getAlternateAlleles().size());
        assertEquals("T", outputVariants.get(1).getAlternateAllele(0).getBaseString());
    }

    // MIVO: variants: 170756001 (C->T),  170756002 (G->T),  170756003 (A->T)
    // MIVO:    4 records contain all 3; 1 record contains 1 and 2; 1 record contains 2 and 3;
    // MIVO:    technically both (1,2) and (2,3) are mnvs, but (1,2,3) is not frequent enough. merge first mnv and output 3rd variant
    @Test
    public void correctOutputFor2PotentialMnvInRegionOf3() {
        final PotentialMNVRegion region = PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(),
                Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12), VARIANTS.get(13)),
                GAP_SIZE);
        final SAMRecord recordWith1and2 = TestUtils.buildSamRecord(170755999, "6M", "CATTAG");
        final SAMRecord recordWith2and3 = TestUtils.buildSamRecord(170755999, "6M", "CACTTG");
        final SAMRecord recordWith1and2and3 = TestUtils.buildSamRecord(170756000, "6M", "ATTTGC");
        final List<SAMRecord> records = Lists.newArrayList(recordWith1and2,
                recordWith1and2and3,
                recordWith1and2and3,
                recordWith1and2and3,
                recordWith1and2and3,
                recordWith2and3);
        final MNVRegionValidator validator = MNVValidator.validateMNVs(records.iterator(), region);
        final List<VariantContext> outputVariants = MNVValidator.outputVariants(validator, MNV_MERGER);
        LOGGER.info(outputVariants);
        assertEquals(2, outputVariants.size());
        assertEquals("CG", outputVariants.get(0).getReference().getBaseString());
        assertEquals(1, outputVariants.get(0).getAlternateAlleles().size());
        assertEquals("TT", outputVariants.get(0).getAlternateAllele(0).getBaseString());
        assertEquals(VARIANTS.get(13), outputVariants.get(1));
    }

    // MIVO: variants: 170756001 (C->T),  170756002 (G->T),  170756003 (A->T),  170756004 (T->C)
    // MIVO:    2 records contain all;
    @Test
    public void correctOutputForMnvOf4InRegionOf4() {
        final PotentialMNVRegion region = PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(),
                Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12), VARIANTS.get(13), VARIANTS.get(14)),
                GAP_SIZE);
        final SAMRecord recordWith1and2and3and4 = TestUtils.buildSamRecord(170756001, "5M", "TTTCC");
        final List<SAMRecord> records = Lists.newArrayList(recordWith1and2and3and4, recordWith1and2and3and4);
        final MNVRegionValidator validator = MNVValidator.validateMNVs(records.iterator(), region);
        final List<VariantContext> outputVariants = MNVValidator.outputVariants(validator, MNV_MERGER);
        assertEquals(1, outputVariants.size());
        assertEquals("CGAT", outputVariants.get(0).getReference().getBaseString());
        assertEquals(1, outputVariants.get(0).getAlternateAlleles().size());
        assertEquals("TTTC", outputVariants.get(0).getAlternateAllele(0).getBaseString());
    }

    // MIVO: variants: 170756001 (C->T),  170756002 (G->T),  170756003 (A->T),  170756004 (T->C)
    // MIVO:    1st record contains 1, 2 and 3 ; 2nd record contains 2, 3 and 4; 4 records contain all
    @Test
    public void correctOutputFor2MnvOf3InRegionOf4() {
        final PotentialMNVRegion region = PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(),
                Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12), VARIANTS.get(13), VARIANTS.get(14)),
                GAP_SIZE);
        final SAMRecord recordWith1and2and3 = TestUtils.buildSamRecord(170755999, "6M", "CATTTT");
        final SAMRecord recordWith2and3and4 = TestUtils.buildSamRecord(170756000, "6M", "ACTTCC");
        final SAMRecord recordWith1and2and3and4 = TestUtils.buildSamRecord(170756001, "5M", "TTTCC");
        final List<SAMRecord> records = Lists.newArrayList(recordWith1and2and3and4,
                recordWith1and2and3and4,
                recordWith1and2and3and4,
                recordWith1and2and3and4,
                recordWith1and2and3,
                recordWith2and3and4);
        final MNVRegionValidator validator = MNVValidator.validateMNVs(records.iterator(), region);
        final List<VariantContext> outputVariants = MNVValidator.outputVariants(validator, MNV_MERGER);
        assertEquals(2, outputVariants.size());
        assertEquals("CGA", outputVariants.get(0).getReference().getBaseString());
        assertEquals(1, outputVariants.get(0).getAlternateAlleles().size());
        assertEquals("TTT", outputVariants.get(0).getAlternateAllele(0).getBaseString());
        assertEquals(VARIANTS.get(14), outputVariants.get(1));
    }

    // MIVO: variants: 170756001 (C->T),  170756002 (G->T),  170756003 (A->T),  170756004 (T->C)
    // MIVO:    1st record contains 1 and 2; 2nd record contains 3 and 4; 3rd record contains all
    @Test
    public void correctOutputFor2MnvOf2InRegionOf4() {
        final PotentialMNVRegion region = PotentialMNVRegion.addVariants(PotentialMNVRegion.empty(),
                Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12), VARIANTS.get(13), VARIANTS.get(14)),
                GAP_SIZE);
        final SAMRecord recordWith1and2 = TestUtils.buildSamRecord(170755999, "6M", "CATTAT");
        final SAMRecord recordWith3and4 = TestUtils.buildSamRecord(170756000, "6M", "ACGTCC");
        final SAMRecord recordWith1and2and3and4 = TestUtils.buildSamRecord(170756001, "5M", "TTTCC");
        final List<SAMRecord> records = Lists.newArrayList(recordWith1and2, recordWith3and4, recordWith1and2and3and4);
        final MNVRegionValidator validator = MNVValidator.validateMNVs(records.iterator(), region);
        final List<VariantContext> outputVariants = MNVValidator.outputVariants(validator, MNV_MERGER);
        assertEquals(2, outputVariants.size());
        assertEquals("CG", outputVariants.get(0).getReference().getBaseString());
        assertEquals(1, outputVariants.get(0).getAlternateAlleles().size());
        assertEquals("TT", outputVariants.get(0).getAlternateAllele(0).getBaseString());
        assertEquals("AT", outputVariants.get(1).getReference().getBaseString());
        assertEquals(1, outputVariants.get(1).getAlternateAlleles().size());
        assertEquals("TC", outputVariants.get(1).getAlternateAllele(0).getBaseString());
    }
}
